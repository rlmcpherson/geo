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

// Polygon represents a sequence of zero or more loops; recall that the
// interior of a loop is defined to be its left-hand side (see Loop).
//
// When the polygon is initialized, the given loops are automatically converted
// into a canonical form consisting of "shells" and "holes".  Shells and holes
// are both oriented CCW, and are nested hierarchically.  The loops are
// reordered to correspond to a preorder traversal of the nesting hierarchy.
//
// Polygons may represent any region of the sphere with a polygonal boundary,
// including the entire sphere (known as the "full" polygon).  The full polygon
// consists of a single full loop (see Loop), whereas the empty polygon has no
// loops at all.
//
// Use FullPolygon() to construct a full polygon. The zero value of Polygon is
// treated as the empty polygon.
//
// Polygons have the following restrictions:
//
//  - Loops may not cross, i.e. the boundary of a loop may not intersect
//    both the interior and exterior of any other loop.
//
//  - Loops may not share edges, i.e. if a loop contains an edge AB, then
//    no other loop may contain AB or BA.
//
//  - Loops may share vertices, however no vertex may appear twice in a
//    single loop (see Loop).
//
//  - No loop may be empty.  The full loop may appear only in the full polygon.
type Polygon struct {
	loops []*Loop
}

// PolygonFromLoops constructs a polygon from the given hierarchically nested
// loops. The polygon interior consists of the points contained by an odd
// number of loops. (Recall that a loop contains the set of points on its
// left-hand side.)
//
// This method figures out the loop nesting hierarchy and assigns every loop a
// depth. Shells have even depths, and holes have odd depths.
//
// NOTE: this function is NOT YET IMPLEMENTED for more than one loop and will
// panic if given a slice of length > 1.
func PolygonFromLoops(loops []*Loop) *Polygon {
	if len(loops) > 1 {
		panic("s2.PolygonFromLoops for multiple loops is not yet implemented")
	}
	return &Polygon{
		loops: loops,
	}
}

// FullPolygon returns a special "full" polygon.
func FullPolygon() *Polygon {
	return &Polygon{
		loops: []*Loop{
			FullLoop(),
		},
	}
}

// IsEmpty reports whether this is the special "empty" polygon (consisting of no loops).
func (p *Polygon) IsEmpty() bool {
	return len(p.loops) == 0
}

// IsFull reports whether this is the special "full" polygon (consisting of a
// single loop that encompasses the entire sphere).
func (p *Polygon) IsFull() bool {
	return len(p.loops) == 1 && p.loops[0].IsFull()
}

// Loops returns the loops in this polygon.
func (p *Polygon) Loops() []*Loop {
	return p.loops
}

// ContainsPoint reports whether this polygon contains the given point
func (p *Polygon) ContainsPoint(pt Point) bool {
	for _, l := range p.loops {
		if l.ContainsPoint(pt) {
			return true
		}
		// TODO: hole check
	}
	return false
}

func (p *Polygon) CapBound() Cap {
	return p.RectBound().CapBound()
}

func (p *Polygon) RectBound() Rect {
	bound := EmptyRect()
	for _, l := range p.loops {
		bound = bound.Union(l.RectBound())
	}
	return bound
}

func (p *Polygon) ContainsCell() bool {
	// TODO: implement for region interface
	return false
}

func (p *Polygon) IntersectsCell() bool {
	// TODO: implement for region interface
	return false
}
