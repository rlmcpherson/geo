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
	"testing"
)

var (
	// The northern hemisphere, defined using two pairs of antipodal points.
	northHemi = LoopFromPoints(parsePoints("0:-180, 0:-90, 0:0, 0:90"))

	// The northern hemisphere, defined using three points 120 degrees apart.
	northHemi3 = LoopFromPoints(parsePoints("0:-180, 0:-60, 0:60"))

	// The southern hemisphere, defined using two pairs of antipodal points.
	southHemi = LoopFromPoints(parsePoints("0:90, 0:0, 0:-90, 0:-180"))

	// The western hemisphere, defined using two pairs of antipodal points.
	westHemi = LoopFromPoints(parsePoints("0:-180, -90:0, 0:0, 90:0"))

	// The eastern hemisphere, defined using two pairs of antipodal points.
	eastHemi = LoopFromPoints(parsePoints("90:0, 0:0, -90:0, 0:-180"))

	// The "near" hemisphere, defined using two pairs of antipodal points.
	nearHemi = LoopFromPoints(parsePoints("0:-90, -90:0, 0:90, 90:0"))

	// The "far" hemisphere, defined using two pairs of antipodal points.
	farHemi = LoopFromPoints(parsePoints("90:0, 0:90, -90:0, 0:-90"))

	// A spiral stripe that slightly over-wraps the equator.
	candyCane = LoopFromPoints(parsePoints("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70"))

	// A small clockwise loop in the northern & eastern hemisperes.
	smallNECW = LoopFromPoints(parsePoints("35:20, 45:20, 40:25"))

	// Loop around the north pole at 80 degrees.
	arctic80 = LoopFromPoints(parsePoints("80:-150, 80:-30, 80:90"))

	// Loop around the south pole at 80 degrees.
	antarctic80 = LoopFromPoints(parsePoints("-80:120, -80:0, -80:-120"))

	// A completely degenerate triangle along the equator that RobustCCW()
	// considers to be CCW.
	lineTriangle = LoopFromPoints(parsePoints("0:1, 0:2, 0:3"))

	// A nearly-degenerate CCW chevron near the equator with very long sides
	// (about 80 degrees).  Its area is less than 1e-640, which is too small
	// to represent in double precision.
	skinnyChevron = LoopFromPoints(parsePoints("0:0, -1e-320:80, 0:1e-320, 1e-320:80"))

	// A diamond-shaped loop around the point 0:180.
	loopA = LoopFromPoints(parsePoints("0:178, -1:180, 0:-179, 1:-180"))

	// Like loopA, but the vertices are at leaf cell centers.
	snappedLoopA = LoopFromPoints([]Point{
		CellIDFromLatLng(parseLatLngs("0:178")[0]).Point(),
		CellIDFromLatLng(parseLatLngs("-1:180")[0]).Point(),
		CellIDFromLatLng(parseLatLngs("0:-179")[0]).Point(),
		CellIDFromLatLng(parseLatLngs("1:-180")[0]).Point(),
	})

	// A different diamond-shaped loop around the point 0:180.
	loopB = LoopFromPoints(parsePoints("0:179, -1:180, 0:-178, 1:-180"))

	// The intersection of A and B.
	aIntersectB = LoopFromPoints(parsePoints("0:179, -1:180, 0:-179, 1:-180"))

	// The union of A and B.
	aUnionB = LoopFromPoints(parsePoints("0:178, -1:180, 0:-178, 1:-180"))

	// A minus B (concave).
	aMinusB = LoopFromPoints(parsePoints("0:178, -1:180, 0:179, 1:-180"))

	// B minus A (concave).
	bMinusA = LoopFromPoints(parsePoints("0:-179, -1:180, 0:-178, 1:-180"))

	// A shape gotten from A by adding a triangle to one edge, and
	// subtracting a triangle from the opposite edge.
	loopC = LoopFromPoints(parsePoints("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180"))

	// A shape gotten from A by adding a triangle to one edge, and
	// adding another triangle to the opposite edge.
	loopD = LoopFromPoints(parsePoints("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180"))

	//   3------------2
	//   |            |               ^
	//   |  7-8  b-c  |               |
	//   |  | |  | |  |      Latitude |
	//   0--6-9--a-d--1               |
	//   |  | |       |               |
	//   |  f-e       |               +----------->
	//   |            |                 Longitude
	//   4------------5
	//
	// Important: It is not okay to skip over collinear vertices when
	// defining these loops (e.g. to define loop E as "0,1,2,3") because S2
	// uses symbolic perturbations to ensure that no three vertices are
	// *ever* considered collinear (e.g., vertices 0, 6, 9 are not
	// collinear).  In other words, it is unpredictable (modulo knowing the
	// details of the symbolic perturbations) whether 0123 contains 06123
	// for example.

	// Loop E:  0,6,9,a,d,1,2,3
	// Loop F:  0,4,5,1,d,a,9,6
	// Loop G:  0,6,7,8,9,a,b,c,d,1,2,3
	// Loop H:  0,6,f,e,9,a,b,c,d,1,2,3
	// Loop I:  7,6,f,e,9,8
	loopE = LoopFromPoints(parsePoints("0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30"))
	loopF = LoopFromPoints(parsePoints("0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34"))
	loopG = LoopFromPoints(parsePoints("0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30"))
	loopH = LoopFromPoints(parsePoints("0:30, 0:34, -10:34, -10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30"))

	loopI = LoopFromPoints(parsePoints("10:34, 0:34, -10:34, -10:36, 0:36, 10:36"))
)

func TestEmptyFullLoops(t *testing.T) {
	emptyLoop := EmptyLoop()

	if !emptyLoop.IsEmpty() {
		t.Errorf("empty loop should be empty")
	}
	if emptyLoop.IsFull() {
		t.Errorf("empty loop should not be full")
	}
	if !emptyLoop.isEmptyOrFull() {
		t.Errorf("empty loop should pass IsEmptyOrFull")
	}

	fullLoop := FullLoop()

	if fullLoop.IsEmpty() {
		t.Errorf("full loop should not be empty")
	}
	if !fullLoop.IsFull() {
		t.Errorf("full loop should be full")
	}
	if !fullLoop.isEmptyOrFull() {
		t.Errorf("full loop should pass IsEmptyOrFull")
	}

	if !empty.RectBound().IsEmpty() {
		t.Errorf("empty loops RectBound should be empty")
	}

	if !full.RectBound().IsFull() {
		t.Errorf("full loops RectBound should be full")
	}

}

func TestLoopContainsPoint(t *testing.T) {

	tests := []struct {
		l    *Loop
		p    Point
		want bool
	}{
		{antarctic80,
			PointFromLatLng(LatLng{southPoleLat, 0}),
			true,
		},
		{southHemi,
			PointFromLatLng(LatLng{southPoleLat, 0}),
			true,
		},
		{antarctic80,
			PointFromLatLng(LatLng{northPoleLat, 0}),
			false,
		},
		{northHemi,
			PointFromLatLng(LatLng{southPoleLat, 0}),
			false,
		},
		{
			arctic80,
			PointFromLatLng(LatLng{northPoleLat, 0}),
			true,
		},
		{
			loopB,
			PointFromLatLng(LatLng{0, 180}),
			true,
		},
	}

	for i, test := range tests {
		// sanity check
		vertices := test.l.Vertices()
		direction := RobustSign(vertices[0], vertices[1], vertices[2])
		if direction == Clockwise {
			t.Errorf("%d: %v not CCW", i, test.l)
		}
		if got := test.l.ContainsPoint(test.p); got != test.want {
			t.Errorf("%d: loop %#v ContainsPoint(%#v), got %v want %v", i, test.l, test.p, got, test.want)
		}

	}

}
