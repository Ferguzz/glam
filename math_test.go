// Copyright (c) 2013 Laurent Moussault. All rights reserved.
// Licensed under a simplified BSD license (see LICENSE file).

package glm

import (
	"math"
	"testing"
)

//------------------------------------------------------------------------------

func TestSqrt(t *testing.T) {
	a := Sqrt(float32(4))
	if a != 2.0 {
		t.Errorf("Wrong result: %v %T", a, a)
	}
	b := Sqrt(float32(3.3))
	if b != 1.8165902 {
		t.Errorf("Wrong result: %v %T", b, b)
	}
}

func BenchmarkMathSqrtFloat64(b *testing.B) {
	a := float64(3.3)
	for i := 0; i < b.N; i++ {
		_ = math.Sqrt(a)
	}
}

func BenchmarkMathSqrtFloat32(b *testing.B) {
	a := float32(3.3)
	for i := 0; i < b.N; i++ {
		_ = float32(math.Sqrt(float64(a)))
	}
}

func BenchmarkGlmSqrt(b *testing.B) {
	a := float32(3.3)
	for i := 0; i < b.N; i++ {
		_ = Sqrt(a)
	}
}

//------------------------------------------------------------------------------

func TestFloor(t *testing.T) {
	a := Floor(float32(3.3))
	if a != 3 {
		t.Errorf("Wrong result for Floor(3.3): %v %T", a, a)
	}
	b := Floor(float32(-3.3))
	if b != -4 {
		t.Errorf("Wrong result for Floor(-3.3): %v %T", b, b)
	}
	c := Floor(float32(3))
	if c != 3 {
		t.Errorf("Wrong result for Floor(3): %v %T", c, c)
	}
	d := Floor(float32(-3))
	if d != -3 {
		t.Errorf("Wrong result for Floor(-3): %v %T", d, d)
	}
}

func BenchmarkMathFloorFloat64(b *testing.B) {
	x := float64(3.3)
	y := float64(-3.3)
	for i := 0; i < b.N; i++ {
		_ = math.Floor(x)
		_ = math.Floor(y)
	}
}

func BenchmarkMathFloorFloat32(b *testing.B) {
	x := float32(3.3)
	y := float32(-3.3)
	for i := 0; i < b.N; i++ {
		_ = float32(math.Floor(float64(x)))
		_ = float32(math.Floor(float64(y)))
	}
}

func BenchmarkGlmFloor(b *testing.B) {
	x := float32(3.3)
	y := float32(-3.3)
	for i := 0; i < b.N; i++ {
		_ = Floor(x)
		_ = Floor(y)
	}
}

//------------------------------------------------------------------------------

func TestFastFloor(t *testing.T) {
	a := FastFloor(float32(3.3))
	if a != 3 {
		t.Errorf("Wrong result for FastFloor(3.3): %v %T", a, a)
	}
	b := FastFloor(float32(-3.3))
	if b != -4 {
		t.Errorf("Wrong result for FastFloor(-3.3): %v %T", b, b)
	}
	c := FastFloor(float32(3))
	if c != 3 {
		t.Errorf("Wrong result for FastFloor(3): %v %T", c, c)
	}
	d := FastFloor(float32(-3))
	if d != -4 {
		t.Errorf("Wrong result for FastFloor(-3): %v %T", d, d)
	}
}

func castFloor(x float32) int32 {
	if x > 0 {
		return int32(x)
	} else {
		return int32(x - 1)
	}
}

func BenchmarkCastFloor(b *testing.B) {
	x := float32(3.3)
	y := float32(-3.3)
	for i := 0; i < b.N; i++ {
		_ = castFloor(x)
		_ = castFloor(y)
	}
}

func BenchmarkGlmFastFloor(b *testing.B) {
	x := float32(3.3)
	y := float32(-3.3)
	for i := 0; i < b.N; i++ {
		_ = FastFloor(x)
		_ = FastFloor(y)
	}
}

//------------------------------------------------------------------------------

func TestRound(t *testing.T) {
	a := Round(float32(3.3))
	if a != 3 {
		t.Errorf("Wrong result for Round(3.3): %v %T", a, a)
	}
	b := Round(float32(-3.3))
	if b != -3 {
		t.Errorf("Wrong result for Round(-3.3): %v %T", b, b)
	}
	c := Round(float32(3))
	if c != 3 {
		t.Errorf("Wrong result for Round(3): %v %T", c, c)
	}
	d := Round(float32(-3))
	if d != -3 {
		t.Errorf("Wrong result for Round(-3): %v %T", d, d)
	}
	e := Round(float32(3.7))
	if e != 4 {
		t.Errorf("Wrong result for Round(3.7): %v %T", e, e)
	}
	f := Round(float32(-3.7))
	if f != -4 {
		t.Errorf("Wrong result for Round(-3.7): %v %T", f, f)
	}
}

func simpleRound(x float32) int32 {
	if x > 0 {
		return int32(x + 0.5)
	} else {
		return int32(x - 0.5)
	}
}

func BenchmarkSimpleRound(b *testing.B) {
	x := float32(3.3)
	y := float32(-3.3)
	for i := 0; i < b.N; i++ {
		_ = simpleRound(x)
		_ = simpleRound(y)
	}
}

func BenchmarkRound(b *testing.B) {
	x := float32(3.3)
	y := float32(-3.3)
	for i := 0; i < b.N; i++ {
		_ = Round(x)
		_ = Round(y)
	}
}

//------------------------------------------------------------------------------
