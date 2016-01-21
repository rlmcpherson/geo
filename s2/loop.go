/*
Copyright 2015 Google Inc. All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

package s2

import (
	"log"

	"github.com/golang/geo/r3"
)

// Loop represents a simple spherical polygon. It consists of a sequence
// of vertices where the first vertex is implicitly connected to the
// last. All loops are defined to have a CCW orientation, i.e. the interior of
// the loop is on the left side of the edges. This implies that a clockwise
// loop enclosing a small area is interpreted to be a CCW loop enclosing a
// very large area.
//
// Loops are not allowed to have any duplicate vertices (whether adjacent or
// not), and non-adjacent edges are not allowed to intersect. Loops must have
// at least 3 vertices (except for the "empty" and "full" loops discussed
// below).
//
// There are two special loops: the "empty" loop contains no points and the
// "full" loop contains all points. These loops do not have any edges, but to
// preserve the invariant that every loop can be represented as a vertex
// chain, they are defined as having exactly one vertex each (see EmptyLoop
// and FullLoop).
type Loop struct {
	vertices []Point

	// originInside keeps a precomputed value whether this loop contains the origin
	// versus computing from the set of vertices every time.
	originInside bool

	// bound is a conservative bound on all points contained by this loop.
	// If l.ContainsPoint(P), then l.bound.ContainsPoint(P).
	bound Rect

	// Since "bound" is not exact, it is possible that a loop A contains
	// another loop B whose bounds are slightly larger. subregionBound
	// has been expanded sufficiently to account for this error, i.e.
	// if A.Contains(B), then A.subregionBound.Contains(B.bound).
	subregionBound Rect
}

// LoopFromPoints constructs a loop from the given points.
func LoopFromPoints(pts []Point) *Loop {
	l := &Loop{
		vertices: pts,
	}

	l.bound = FullRect()
	l.initOriginAndBound()
	return l
}

// EmptyLoop returns a special "empty" loop.
func EmptyLoop() *Loop {
	return LoopFromPoints([]Point{Point{r3.Vector{X: 0, Y: 0, Z: 1}}})
}

// FullLoop returns a special "full" loop.
func FullLoop() *Loop {
	return LoopFromPoints([]Point{Point{r3.Vector{X: 0, Y: 0, Z: -1}}})
}

// initOriginAndBound sets the origin containment for the given point and then calls
// the initialization for the bounds objects and the internal index.
func (l *Loop) initOriginAndBound() {
	if len(l.vertices) < 3 {
		// Check for the special "empty" and "full" loops (which have one vertex).
		if !l.isEmptyOrFull() {
			l.originInside = false
			return
		}

		// This is the special empty or full loop, so the origin depends on if
		// the vertex is in the southern hemisphere or not.
		l.originInside = l.vertices[0].Z < 0
	} else {
		if !l.bound.ContainsPoint(l.vertices[1]) {
			panic("bound must contain vertex 1")
		}

		v1Inside := OrderedCCW(Point{l.vertices[1].Ortho()}, l.vertices[0], l.vertices[2], l.vertices[1])
		if v1Inside != l.ContainsPoint(l.vertices[1]) {
			l.originInside = true
		}
	}

	// We *must* call initBound() before initIndex(), because initBound() calls
	// Contains(s2.Point), and Contains(s2.Point) does a bounds check whenever the
	// index is not fresh (i.e., the loop has been added to the index but the
	// index has not been updated yet).
	l.initBound()

	// TODO(roberts): Depends on s2shapeindex being implemented.
	// l.initIndex()
}

// initBound sets up the approximate bounding Rects for this loop.
//
// TODO(roberts): Full implementation depends on s2edgeutil and s2rectbounder.
func (l *Loop) initBound() {
	// Check for the special "empty" and "full" loops.
	if l.isEmptyOrFull() {
		if l.IsEmpty() {
			l.bound = Rect{}
		} else {
			l.bound = FullRect()
		}
		l.subregionBound = l.bound
		return
	}

	// TODO: check init with 0 vertices
	bound := EmptyRect()
	for _, p := range l.vertices {
		bound = bound.AddPoint(LatLngFromPoint(p))
	}
	l.bound, l.subregionBound = bound, bound

	// TODO(roberts): As other s2 pieces get ported, complete this method.
	return
}

// ContainsOrigin reports true if this loop contains s2.Origin().
func (l Loop) ContainsOrigin() bool {
	return l.originInside
}

// IsEmpty reports true if this is the special "empty" loop that contains no points.
func (l Loop) IsEmpty() bool {
	return l.isEmptyOrFull() && !l.ContainsOrigin()
}

// IsFull reports true if this is the special "full" loop that contains all points.
func (l Loop) IsFull() bool {
	return l.isEmptyOrFull() && l.ContainsOrigin()
}

// isEmptyOrFull reports true if this loop is either the "empty" or "full" special loops.
func (l Loop) isEmptyOrFull() bool {
	return len(l.vertices) == 1
}

// RectBound returns a tight bounding rectangle. If the loop contains the point,
// the bound also contains it.
func (l Loop) RectBound() Rect {
	return l.bound
}

// Vertices returns the vertices in the loop.
func (l Loop) Vertices() []Point {
	return l.vertices
}

// ContainsPoint reports whether this loop contains the given point
func (l Loop) ContainsPoint(p Point) bool {
	log.Printf("bound hi: %#v, lo: %#v", l.RectBound().Hi(), l.RectBound().Lo())
	if !l.RectBound().ContainsPoint(p) {
		return false
	}
	log.Println("point in rect")
	contains := l.originInside
	v0, origin := l.vertices[0], OriginPoint()
	for _, v := range l.vertices[1:] {
		contains = contains != EdgeOrVertexCrossing(origin, p, v0, v)
	}
	return contains
}

// BUG(): The major differences from the C++ version is pretty much everything.
