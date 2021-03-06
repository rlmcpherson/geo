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
	"math"
	"testing"
)

func TestCol(t *testing.T) {
	tests := []struct {
		have *matrix3x3
		col  int
		want Point
	}{
		{&matrix3x3{}, 0, OriginPoint()},
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			0,
			PointFromCoords(1, 4, 7),
		},
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			2,
			PointFromCoords(3, 6, 9),
		},
	}

	for _, test := range tests {
		if got := test.have.col(test.col); !got.ApproxEqual(test.want) {
			t.Errorf("%v.col(%d) = %v, want %v", test.have, test.col, got, test.want)
		}
	}
}

func TestRow(t *testing.T) {
	tests := []struct {
		have *matrix3x3
		row  int
		want Point
	}{
		{&matrix3x3{}, 0, OriginPoint()},
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			0,
			PointFromCoords(1, 2, 3),
		},
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			2,
			PointFromCoords(7, 8, 9),
		},
	}

	for _, test := range tests {
		if got := test.have.row(test.row); !got.ApproxEqual(test.want) {
			t.Errorf("%v.row(%d) = %v, want %v", test.have, test.row, got, test.want)
		}
	}
}

func TestSetCol(t *testing.T) {
	tests := []struct {
		have  *matrix3x3
		col   int
		point Point
		want  *matrix3x3
	}{
		{
			&matrix3x3{},
			0,
			PointFromCoords(1, 1, 0),
			&matrix3x3{
				{math.Sqrt(2) / 2, 0, 0},
				{math.Sqrt(2) / 2, 0, 0},
				{0, 0, 0},
			},
		},
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			2,
			PointFromCoords(1, 1, 0),
			&matrix3x3{
				{1, 2, math.Sqrt(2) / 2},
				{4, 5, math.Sqrt(2) / 2},
				{7, 8, 0},
			},
		},
	}

	for _, test := range tests {
		if got := test.have.setCol(test.col, test.point); !matricesApproxEqual(got, test.want) {
			t.Errorf("%v.setCol(%d, %v) = %v, want %v", test.have, test.col, test.point, got, test.want)
		}
	}
}

func TestSetRow(t *testing.T) {
	tests := []struct {
		have  *matrix3x3
		row   int
		point Point
		want  *matrix3x3
	}{
		{
			&matrix3x3{},
			0,
			PointFromCoords(1, 1, 0),
			&matrix3x3{
				{math.Sqrt(2) / 2, math.Sqrt(2) / 2, 0},
				{0, 0, 0},
				{0, 0, 0},
			},
		},
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			2,
			PointFromCoords(1, 1, 0),
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{math.Sqrt(2) / 2, math.Sqrt(2) / 2, 0},
			},
		},
	}
	for _, test := range tests {
		if got := test.have.setRow(test.row, test.point); !matricesApproxEqual(got, test.want) {
			t.Errorf("%v.setRow(%d, %v) = %v, want %v", test.have, test.row, test.point, got, test.want)
		}
	}
}

func TestScale(t *testing.T) {
	tests := []struct {
		have  *matrix3x3
		scale float64
		want  *matrix3x3
	}{
		{
			&matrix3x3{},
			0,
			&matrix3x3{},
		},
		{
			&matrix3x3{
				{1, 1, 1},
				{1, 1, 1},
				{1, 1, 1},
			},
			0,
			&matrix3x3{},
		},
		{
			&matrix3x3{
				{1, 1, 1},
				{1, 1, 1},
				{1, 1, 1},
			},
			1,
			&matrix3x3{
				{1, 1, 1},
				{1, 1, 1},
				{1, 1, 1},
			},
		},
		{
			&matrix3x3{
				{1, 1, 1},
				{1, 1, 1},
				{1, 1, 1},
			},
			5,
			&matrix3x3{
				{5, 5, 5},
				{5, 5, 5},
				{5, 5, 5},
			},
		},
		{
			&matrix3x3{
				{-2, 2, -3},
				{-1, 1, 3},
				{2, 0, -1},
			},
			2.75,
			&matrix3x3{
				{-5.5, 5.5, -8.25},
				{-2.75, 2.75, 8.25},
				{5.5, 0, -2.75},
			},
		},
	}

	for _, test := range tests {
		if got := test.have.scale(test.scale); !matricesApproxEqual(got, test.want) {
			t.Errorf("%v.scale(%f) = %v, want %v", test.have, test.scale, got, test.want)
		}
	}
}

func TestMul(t *testing.T) {
	tests := []struct {
		have  *matrix3x3
		point Point
		want  Point
	}{
		{&matrix3x3{}, Point{}, Point{}},
		{
			&matrix3x3{
				{1, 1, 1},
				{1, 1, 1},
				{1, 1, 1},
			},
			Point{},
			Point{},
		},
		{
			// Identity times something gives back the something
			&matrix3x3{
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1},
			},
			Point{},
			Point{},
		},
		{
			// Identity times something gives back the something
			&matrix3x3{
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1},
			},
			PointFromCoords(1, 2, 3),
			PointFromCoords(1, 2, 3),
		},
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			PointFromCoords(1, 1, 1),
			PointFromCoords(6, 15, 24),
		},
	}
	for _, test := range tests {
		if got := test.have.mul(test.point); !got.ApproxEqual(test.want) {
			t.Errorf("%v.mul(%v) = %v, want %v", test.have, test.point, got, test.want)
		}
	}
}

func TestDet(t *testing.T) {
	tests := []struct {
		have *matrix3x3
		want float64
	}{
		{
			&matrix3x3{},
			0,
		},
		{
			// Matrix of all the same values has det of 0.
			&matrix3x3{
				{1, 1, 1},
				{1, 1, 1},
				{1, 1, 1},
			},
			0,
		},
		{
			// Identity matrix has det of 1.
			&matrix3x3{
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1},
			},
			1,
		},
		{
			&matrix3x3{
				{-2, 2, -3},
				{-1, 1, 3},
				{2, 0, -1},
			},
			18,
		},
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			0,
		},
		{
			&matrix3x3{
				{9, 8, 7},
				{6, 5, 4},
				{3, 2, 1},
			},
			0,
		},
		{
			&matrix3x3{
				{1.74, math.E, 42},
				{math.Pi, math.Sqrt2, math.Ln10},
				{3, math.SqrtPhi, 9.8976},
			},
			-56.838525224123096,
		},
	}

	for _, test := range tests {
		if got := test.have.det(); !float64Eq(got, test.want) {
			t.Errorf("%v.det() = %v, want %v", test.have, got, test.want)
		}
	}
}

func TestTranspose(t *testing.T) {
	tests := []struct {
		have *matrix3x3
		want *matrix3x3
	}{
		{&matrix3x3{}, &matrix3x3{}},
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			&matrix3x3{
				{1, 4, 7},
				{2, 5, 8},
				{3, 6, 9},
			},
		},
		{
			&matrix3x3{
				{1, 0, 0},
				{0, 2, 0},
				{0, 0, 3},
			},
			&matrix3x3{
				{1, 0, 0},
				{0, 2, 0},
				{0, 0, 3},
			},
		},
		{
			&matrix3x3{
				{1, 2, 3},
				{0, 4, 5},
				{0, 0, 6},
			},
			&matrix3x3{
				{1, 0, 0},
				{2, 4, 0},
				{3, 5, 6},
			},
		},
		{
			&matrix3x3{
				{1, 1, 1},
				{0, 0, 0},
				{0, 0, 0},
			},
			&matrix3x3{
				{1, 0, 0},
				{1, 0, 0},
				{1, 0, 0},
			},
		},
	}

	for _, test := range tests {
		if got := test.have.transpose().transpose(); !matricesApproxEqual(got, test.have) {
			t.Errorf("%v.transpose().transpose() = %v, want %v", test.have, got, test.have)
		}

		if got := test.have.transpose(); !matricesApproxEqual(got, test.want) {
			t.Errorf("%v.transpose() = %v, want %v", test.have, got, test.want)
		}

	}
}

func TestString(t *testing.T) {
	tests := []struct {
		have *matrix3x3
		want string
	}{
		{
			&matrix3x3{
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
			},
			`[ 1.0000 2.0000 3.0000 ] [ 4.0000 5.0000 6.0000 ] [ 7.0000 8.0000 9.0000 ]`,
		},
		{
			&matrix3x3{
				{1, 4, 7},
				{2, 5, 8},
				{3, 6, 9},
			},
			`[ 1.0000 4.0000 7.0000 ] [ 2.0000 5.0000 8.0000 ] [ 3.0000 6.0000 9.0000 ]`,
		},
	}

	for _, test := range tests {
		if got := test.have.String(); got != test.want {
			t.Errorf("%v.String() = %v, want %v", test.have, got, test.want)
		}
	}
}
