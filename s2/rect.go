/*
Copyright 2014 Google Inc. All rights reserved.

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
	"fmt"
	"math"

	"github.com/golang/geo/r1"
	"github.com/golang/geo/s1"
)

// Rect represents a closed latitude-longitude rectangle.
type Rect struct {
	Lat r1.Interval
	Lng s1.Interval
}

var (
	validRectLatRange = r1.Interval{-math.Pi / 2, math.Pi / 2}
	validRectLngRange = s1.FullInterval()
)

// EmptyRect returns the empty rectangle.
func EmptyRect() Rect { return Rect{r1.EmptyInterval(), s1.EmptyInterval()} }

// FullRect returns the full rectangle.
func FullRect() Rect { return Rect{validRectLatRange, validRectLngRange} }

// RectFromLatLng constructs a rectangle containing a single point p.
func RectFromLatLng(p LatLng) Rect {
	return Rect{
		Lat: r1.Interval{p.Lat.Radians(), p.Lat.Radians()},
		Lng: s1.Interval{p.Lng.Radians(), p.Lng.Radians()},
	}
}

// RectFromCenterSize constructs a rectangle with the given size and center.
// center needs to be normalized, but size does not. The latitude
// interval of the result is clamped to [-90,90] degrees, and the longitude
// interval of the result is FullRect() if and only if the longitude size is
// 360 degrees or more.
//
// Examples of clamping (in degrees):
//   center=(80,170),  size=(40,60)   -> lat=[60,90],   lng=[140,-160]
//   center=(10,40),   size=(210,400) -> lat=[-90,90],  lng=[-180,180]
//   center=(-90,180), size=(20,50)   -> lat=[-90,-80], lng=[155,-155]
func RectFromCenterSize(center, size LatLng) Rect {
	half := LatLng{size.Lat / 2, size.Lng / 2}
	return RectFromLatLng(center).expanded(half)
}

// IsValid returns true iff the rectangle is valid.
// This requires Lat ⊆ [-π/2,π/2] and Lng ⊆ [-π,π], and Lat = ∅ iff Lng = ∅
func (r Rect) IsValid() bool {
	return math.Abs(r.Lat.Lo) <= math.Pi/2 &&
		math.Abs(r.Lat.Hi) <= math.Pi/2 &&
		r.Lng.IsValid() &&
		r.Lat.IsEmpty() == r.Lng.IsEmpty()
}

// IsEmpty reports whether the rectangle is empty.
func (r Rect) IsEmpty() bool { return r.Lat.IsEmpty() }

// IsFull reports whether the rectangle is full.
func (r Rect) IsFull() bool { return r.Lat.Equal(validRectLatRange) && r.Lng.IsFull() }

// IsPoint reports whether the rectangle is a single point.
func (r Rect) IsPoint() bool { return r.Lat.Lo == r.Lat.Hi && r.Lng.Lo == r.Lng.Hi }

// Lo returns one corner of the rectangle.
func (r Rect) Lo() LatLng {
	return LatLng{s1.Angle(r.Lat.Lo) * s1.Radian, s1.Angle(r.Lng.Lo) * s1.Radian}
}

// Hi returns the other corner of the rectangle.
func (r Rect) Hi() LatLng {
	return LatLng{s1.Angle(r.Lat.Hi) * s1.Radian, s1.Angle(r.Lng.Hi) * s1.Radian}
}

// Center returns the center of the rectangle.
func (r Rect) Center() LatLng {
	return LatLng{s1.Angle(r.Lat.Center()) * s1.Radian, s1.Angle(r.Lng.Center()) * s1.Radian}
}

// Size returns the size of the Rect.
func (r Rect) Size() LatLng {
	return LatLng{s1.Angle(r.Lat.Length()) * s1.Radian, s1.Angle(r.Lng.Length()) * s1.Radian}
}

// Area returns the surface area of the Rect.
func (r Rect) Area() float64 {
	if r.IsEmpty() {
		return 0
	}
	capDiff := math.Abs(math.Sin(r.Lat.Hi) - math.Sin(r.Lat.Lo))
	return r.Lng.Length() * capDiff
}

// ContainsLatLng reports whether the given LatLng is within the Rect.
func (r Rect) ContainsLatLng(ll LatLng) bool {
	if !ll.IsValid() {
		return false
	}
	return r.Lat.Contains(ll.Lat.Radians()) && r.Lng.Contains(ll.Lng.Radians())
}

// AddPoint increases the size of the rectangle to include the given point.
func (r Rect) AddPoint(ll LatLng) Rect {
	if !ll.IsValid() {
		return r
	}
	return Rect{
		Lat: r.Lat.AddPoint(ll.Lat.Radians()),
		Lng: r.Lng.AddPoint(ll.Lng.Radians()),
	}
}

// expanded returns a rectangle that has been expanded by margin.Lat on each side
// in the latitude direction, and by margin.Lng on each side in the longitude
// direction. If either margin is negative, then it shrinks the rectangle on
// the corresponding sides instead. The resulting rectangle may be empty.
//
// The latitude-longitude space has the topology of a cylinder. Longitudes
// "wrap around" at +/-180 degrees, while latitudes are clamped to range [-90, 90].
// This means that any expansion (positive or negative) of the full longitude range
// remains full (since the "rectangle" is actually a continuous band around the
// cylinder), while expansion of the full latitude range remains full only if the
// margin is positive.
//
// If either the latitude or longitude interval becomes empty after
// expansion by a negative margin, the result is empty.
//
// Note that if an expanded rectangle contains a pole, it may not contain
// all possible lat/lng representations of that pole, e.g., both points [π/2,0]
// and [π/2,1] represent the same pole, but they might not be contained by the
// same Rect.
//
// If you are trying to grow a rectangle by a certain distance on the
// sphere (e.g. 5km), refer to the ExpandedByDistance() C++ method implementation
// instead.
func (r Rect) expanded(margin LatLng) Rect {
	lat := r.Lat.Expanded(margin.Lat.Radians())
	lng := r.Lng.Expanded(margin.Lng.Radians())

	if lat.IsEmpty() || lng.IsEmpty() {
		return EmptyRect()
	}

	return Rect{
		Lat: lat.Intersection(validRectLatRange),
		Lng: lng,
	}
}

func (r Rect) String() string { return fmt.Sprintf("[Lo%v, Hi%v]", r.Lo(), r.Hi()) }

// PolarClosure returns the rectangle unmodified if it does not include either pole.
// If it includes either pole, PolarClosure returns an expansion of the rectangle along
// the longitudinal range to include all possible representations of the contained poles.
func (r Rect) PolarClosure() Rect {
	if r.Lat.Lo == -math.Pi/2 || r.Lat.Hi == math.Pi/2 {
		return Rect{r.Lat, s1.FullInterval()}
	}
	return r
}

// Union returns the smallest Rect containing the union of this rectangle and the given rectangle.
func (r Rect) Union(other Rect) Rect {
	return Rect{
		Lat: r.Lat.Union(other.Lat),
		Lng: r.Lng.Union(other.Lng),
	}
}

// Intersection returns the smallest rectangle containing the intersection of
// this rectangle and the given rectangle. Note that the region of intersection
// may consist of two disjoint rectangles, in which case a single rectangle
// spanning both of them is returned.
func (r Rect) Intersection(other Rect) Rect {
	lat := r.Lat.Intersection(other.Lat)
	lng := r.Lng.Intersection(other.Lng)

	if lat.IsEmpty() || lng.IsEmpty() {
		return EmptyRect()
	}
	return Rect{lat, lng}
}

// Intersects reports whether this rectangle and the other have any points in common.
func (r Rect) Intersects(other Rect) bool {
	return r.Lat.Intersects(other.Lat) && r.Lng.Intersects(other.Lng)
}

// CapBound returns a cap that countains Rect.
func (r Rect) CapBound() Cap {
	// We consider two possible bounding caps, one whose axis passes
	// through the center of the lat-long rectangle and one whose axis
	// is the north or south pole.  We return the smaller of the two caps.

	if r.IsEmpty() {
		return EmptyCap()
	}

	var poleZ, poleAngle float64
	if r.Lat.Hi+r.Lat.Lo < 0 {
		// South pole axis yields smaller cap.
		poleZ = -1
		poleAngle = math.Pi/2 + r.Lat.Hi
	} else {
		poleZ = 1
		poleAngle = math.Pi/2 - r.Lat.Lo
	}
	poleCap := CapFromCenterAngle(PointFromCoords(0, 0, poleZ), s1.Angle(poleAngle)*s1.Radian)

	// For bounding rectangles that span 180 degrees or less in longitude, the
	// maximum cap size is achieved at one of the rectangle vertices.  For
	// rectangles that are larger than 180 degrees, we punt and always return a
	// bounding cap centered at one of the two poles.
	if math.Remainder(r.Lng.Hi-r.Lng.Lo, 2*math.Pi) >= 0 && r.Lng.Hi-r.Lng.Lo < 2*math.Pi {
		midCap := CapFromPoint(PointFromLatLng(r.Center())).AddPoint(PointFromLatLng(r.Lo())).AddPoint(PointFromLatLng(r.Hi()))
		if midCap.Height() < poleCap.Height() {
			return midCap
		}
	}
	return poleCap
}

// RectBound returns itself.
func (r Rect) RectBound() Rect {
	return r
}

// Contains reports whether this Rect contains the other Rect.
func (r Rect) Contains(other Rect) bool {
	return r.Lat.ContainsInterval(other.Lat) && r.Lng.ContainsInterval(other.Lng)
}

// ContainsCell reports whether the given point is contained by this Rect.
func (r Rect) ContainsCell(c Cell) bool {
	// A latitude-longitude rectangle contains a cell if and only if it contains
	// the cell's bounding rectangle. This test is exact from a mathematical
	// point of view, assuming that the bounds returned by Cell.RectBound()
	// are tight. However, note that there can be a loss of precision when
	// converting between representations -- for example, if an s2.Cell is
	// converted to a polygon, the polygon's bounding rectangle may not contain
	// the cell's bounding rectangle. This has some slightly unexpected side
	// effects; for instance, if one creates an s2.Polygon from an s2.Cell, the
	// polygon will contain the cell, but the polygon's bounding box will not.
	return r.Contains(c.RectBound())
}

func (r Rect) ContainsPoint(p Point) bool {
	return r.ContainsLatLng(LatLngFromPoint(p))
}

// BUG(dsymonds): The major differences from the C++ version are:
//   - almost everything
