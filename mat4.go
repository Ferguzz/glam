// Copyright (c) 2013 Laurent Moussault. All rights reserved.
// Licensed under a simplified BSD license (see LICENSE file).

package glam

import "github.com/Ferguzz/glam/math"

//------------------------------------------------------------------------------

// `Mat4` is a single-precision matrix with 4 columns and 4 rows.
//
// Note: matrices are stored in column-major order, so when writing literals
// remember to use the transpose.
type Mat4 [16]float32

//------------------------------------------------------------------------------

// `NewMat4` allocates and returns a new matrix. The elements are stored in
// alphabetical order (column-major order).
//
// See also `MakeMat4` and `SetTo`.
func NewMat4(
	a, e, i, m,
	b, f, j, n,
	c, g, k, o,
	d, h, l, p float32,
) *Mat4 {
	return &Mat4{
		a, b, c, d,
		e, f, g, h,
		i, j, k, l,
		m, n, o, p,
	}
}

// `MakeMat4` returns a matrix. The elements are stored in
// alphabetical order (column-major order).
//
// See also `NewMat4` and `SetTo`.
func MakeMat4(
	a, e, i, m,
	b, f, j, n,
	c, g, k, o,
	d, h, l, p float32,
) Mat4 {
	return Mat4{
		a, b, c, d,
		e, f, g, h,
		i, j, k, l,
		m, n, o, p,
	}
}

// `SetTo` initializes `matrix`. The elements are stored in
// alphabetical order (column-major order).
//
// See also `NewMat4` and `SetTo`.
func (matrix *Mat4) SetTo(
	a, e, i, m,
	b, f, j, n,
	c, g, k, o,
	d, h, l, p float32,
) {
	matrix[0] = a
	matrix[1] = b
	matrix[2] = c
	matrix[3] = d

	matrix[4] = e
	matrix[5] = f
	matrix[6] = g
	matrix[7] = h

	matrix[8] = i
	matrix[9] = j
	matrix[10] = k
	matrix[11] = l

	matrix[12] = m
	matrix[13] = n
	matrix[14] = o
	matrix[15] = p
}

//------------------------------------------------------------------------------

// `Zeros` returns a 4x4 zeroed matrix.
func Zeros() Mat4 {
	return MakeMat4(
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0,
	)
}

// `Identity` returns a 4x4 Identity matrix.
func Identity() Mat4 {
	return MakeMat4(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1,
	)
}

//------------------------------------------------------------------------------

// // `At` returns the element at '(row, column)`.
// func (m Mat4) At(row, column int) float32 {
// 	return m[column][row]
// }

// // `Set` sets the element at `(row, column)` to `value`.
// func (m *Mat4) Set(row, column int, value float32) {
// 	m[column][row] = value
// }

//------------------------------------------------------------------------------

// `Perspective` returns a perspective projection matrix.
//
// See also `SetToPerspective`, `PerspectiveFrustum` and `SetToPerspectiveFrustum`.
func Perspective(fieldOfView, aspectRatio, near, far float32) Mat4 {
	f := float32(1.0) / math.Tan(fieldOfView/float32(2.0))

	return Mat4{
		f / aspectRatio, 0, 0, 0,
		0, f, 0, 0,
		0, 0, (far + near) / (near - far), -1,
		0, 0, (2 * far * near) / (near - far), 0,
	}
}

// `SetToPerspective` sets `m` to a perspective projection matrix.
//
// See also `Perspective`, `PerspectiveFrustum` and `SetToPerspectiveFrustum`.
func (m *Mat4) SetToPerspective(fieldOfView, aspectRatio, near, far float32) {
	f := float32(1.0) / math.Tan(fieldOfView/float32(2.0))

	m[0] = f / aspectRatio
	m[1] = 0
	m[2] = 0
	m[3] = 0

	m[4] = 0
	m[5] = f
	m[6] = 0
	m[7] = 0

	m[8] = 0
	m[9] = 0
	m[10] = (far + near) / (near - far)
	m[11] = -1

	m[12] = 0
	m[13] = 0
	m[14] = (2 * far * near) / (near - far)
	m[15] = 0
}

//------------------------------------------------------------------------------

// `PerspectiveFrustum` returns a perspective projection matrix.
//
// See also `SetToPerspectiveFrustum`, `Perspective` and `SetToPerspective`.
func PerspectiveFrustum(left, right, bottom, top, near, far float32) Mat4 {
	return Mat4{
		(2 * near) / (right - left), 0, 0, 0,
		0, (2 * near) / (top - bottom), 0, 0,
		(right + left) / (right - left), (top + bottom) / (top - bottom), -(far + near) / (far - near), -1,
		0, 0, -(2 * far * near) / (far - near), 0,
	}
}

// `SetToPerspectiveFrustum` sets `m` to a perspective projection matrix.
//
// See also `PerspectiveFrustum`, `Perspective` and `SetToPerspective`.
func (m *Mat4) SetToPerspectiveFrustum(left, right, bottom, top, near, far float32) {
	m[0] = (2 * near) / (right - left)
	m[1] = 0
	m[2] = 0
	m[3] = 0

	m[4] = 0
	m[5] = (2 * near) / (top - bottom)
	m[6] = 0
	m[7] = 0

	m[8] = (right + left) / (right - left)
	m[9] = (top + bottom) / (top - bottom)
	m[10] = -(far + near) / (far - near)
	m[11] = -1

	m[12] = 0
	m[13] = 0
	m[14] = -(2 * far * near) / (far - near)
	m[15] = 0
}

//------------------------------------------------------------------------------

// `Orthographic` returns an orthographic (parallel) projection matrix.
// `zoom` is the height of the projection plane.
//
// See also `SetToOrthographic`, `OrthographicFrustum` and `SetToOrthographicFrustum`.
func Orthographic(zoom, aspectRatio, near, far float32) Mat4 {
	top := zoom / 2
	right := top * aspectRatio
	return Mat4{
		1 / right, 0, 0, 0,
		0, 1 / top, 0, 0,
		0, 0, -2 / (far - near), 0,
		0, 0, -(far + near) / (far - near), 1,
	}
}

// `SetToOrthographic` sets `m` to an orthographic (parallel) projection matrix.
// `zoom` is the height of the projection plane.
//
// See also `Orthographic`, `OrthographicFrustum` and `SetToOrthographicFrustum`.
func (m *Mat4) SetToOrthographic(zoom, aspectRatio, near, far float32) {
	top := zoom / 2
	right := top * aspectRatio

	m[0] = 1 / right
	m[1] = 0
	m[2] = 0
	m[3] = 0

	m[4] = 0
	m[5] = 1 / top
	m[6] = 0
	m[7] = 0

	m[8] = 0
	m[9] = 0
	m[10] = -2 / (far - near)
	m[11] = 0

	m[12] = 0
	m[13] = 0
	m[14] = -(far + near) / (far - near)
	m[15] = 1
}

//------------------------------------------------------------------------------

// `OrthographicFrustum` returns an orthographic (parallel) projection matrix.
//
// See also `SetToOrthographicFrustum`, `Orthographic` and `SetToOrthographic`.
func OrthographicFrustum(left, right, bottom, top, near, far float32) Mat4 {
	return Mat4{
		2 / (right - left), 0, 0, 0,
		0, 2 / (top - bottom), 0, 0,
		0, 0, -2 / (far - near), 0,
		-(right + left) / (right - left), -(top + bottom) / (top - bottom), -(far + near) / (far - near), 1,
	}
}

// `SetToOrthographicFrustum` returns an orthographic (parallel) projection matrix.
//
// See also `OrthographicFrustum`, `Orthographic` and `SetToOrthographic`.
func (m *Mat4) SetToOrthographicFrustum(left, right, bottom, top, near, far float32) {
	m[0] = 2 / (right - left)
	m[1] = 0
	m[2] = 0
	m[3] = 0

	m[4] = 0
	m[5] = 2 / (top - bottom)
	m[6] = 0
	m[7] = 0

	m[8] = 0
	m[9] = 0
	m[10] = -2 / (far - near)
	m[11] = 0

	m[12] = -(right + left) / (right - left)
	m[13] = -(top + bottom) / (top - bottom)
	m[14] = -(far + near) / (far - near)
	m[15] = 1
}

//------------------------------------------------------------------------------

// `Translation` returns a translation matrix.
//
// See also `SetToTranslation`.
func Translation(t Vec3) Mat4 {
	return Mat4{
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		t.X, t.Y, t.Z, 1,
	}
}

// `SetToTranslation` sets `m` to a translation matrix.
//
// See also `Translation`.
func (m *Mat4) SetToTranslation(t Vec3) {
	m[0] = 1
	m[1] = 0
	m[2] = 0
	m[3] = 0

	m[4] = 0
	m[5] = 1
	m[6] = 0
	m[7] = 0

	m[8] = 0
	m[9] = 0
	m[10] = 1
	m[11] = 0

	m[12] = t.X
	m[13] = t.Y
	m[14] = t.Z
	m[15] = 1
}

//------------------------------------------------------------------------------

// `Rotation` returns a rotation matrix.
//
// See also `SetToRotation`.
func Rotation(angle float32, axis Vec3) Mat4 {
	c := math.Cos(angle)
	s := math.Sin(angle)

	return Mat4{
		c + axis.X*axis.X*(1-c), -axis.Z*s + axis.X*axis.Y*(1-c), axis.Y*s + axis.X*axis.Z*(1-c), 0,
		axis.Z*s + axis.Y*axis.X*(1-c), c + axis.Y*axis.Y*(1-c), -axis.X*s + axis.Y*axis.Z*(1-c), 0,
		-axis.Y*s + axis.Z*axis.X*(1-c), axis.X*s + axis.Z*axis.Y*(1-c), c + axis.Z*axis.Z*(1-c), 0,
		0, 0, 0, 1,
	}
}

// `SetToRotation` sets `m` to a rotation matrix.
//
// See also `Rotation`.
func (m *Mat4) Rotation(angle float32, axis Vec3) {
	c := math.Cos(angle)
	s := math.Sin(angle)

	m[0] = c + axis.X*axis.X*(1-c)
	m[1] = -axis.Z*s + axis.X*axis.Y*(1-c)
	m[2] = axis.Y*s + axis.X*axis.Z*(1-c)
	m[3] = 0

	m[4] = axis.Z*s + axis.Y*axis.X*(1-c)
	m[5] = c + axis.Y*axis.Y*(1-c)
	m[6] = -axis.X*s + axis.Y*axis.Z*(1-c)
	m[7] = 0

	m[8] = -axis.Y*s + axis.Z*axis.X*(1-c)
	m[9] = axis.X*s + axis.Z*axis.Y*(1-c)
	m[10] = c + axis.Z*axis.Z*(1-c)
	m[11] = 0

	m[12] = 0
	m[13] = 0
	m[14] = 0
	m[15] = 1
}

//------------------------------------------------------------------------------

// `Times` returns the matrix product of `m` and `o`.
//
// See also `Multiply`, `TimesVec` and `MultiplyVec`.
func (m *Mat4) Times(o *Mat4) Mat4 {
	return Mat4{
		m[0]*o[0] + m[1]*o[4] + m[2]*o[8] + m[3]*o[12],
		m[0]*o[1] + m[1]*o[5] + m[2]*o[9] + m[3]*o[13],
		m[0]*o[2] + m[1]*o[6] + m[2]*o[10] + m[3]*o[14],
		m[0]*o[3] + m[1]*o[7] + m[2]*o[11] + m[3]*o[15],

		m[4]*o[0] + m[5]*o[4] + m[6]*o[8] + m[7]*o[12],
		m[4]*o[1] + m[5]*o[5] + m[6]*o[9] + m[7]*o[13],
		m[4]*o[2] + m[5]*o[6] + m[6]*o[10] + m[7]*o[14],
		m[4]*o[3] + m[5]*o[7] + m[6]*o[11] + m[7]*o[15],

		m[8]*o[0] + m[9]*o[4] + m[10]*o[8] + m[11]*o[12],
		m[8]*o[1] + m[9]*o[5] + m[10]*o[9] + m[11]*o[13],
		m[8]*o[2] + m[9]*o[6] + m[10]*o[10] + m[11]*o[14],
		m[8]*o[3] + m[9]*o[7] + m[10]*o[11] + m[11]*o[15],

		m[12]*o[0] + m[13]*o[4] + m[14]*o[8] + m[15]*o[12],
		m[12]*o[1] + m[13]*o[5] + m[14]*o[9] + m[15]*o[13],
		m[12]*o[2] + m[13]*o[6] + m[14]*o[10] + m[15]*o[14],
		m[12]*o[3] + m[13]*o[7] + m[14]*o[11] + m[15]*o[15],
	}
}

// `Multiply` sets `r` to the matrix product of `m` and `o`.
//
// `r` must not be `m` or `o`.
//
// See also `Multiply`, `TimesVec` and `MultiplyVec`.
func (r *Mat4) Multiply(m, o *Mat4) {
	r[0] = m[0]*o[0] + m[1]*o[4] + m[2]*o[8] + m[3]*o[12]
	r[1] = m[0]*o[1] + m[1]*o[5] + m[2]*o[9] + m[3]*o[13]
	r[2] = m[0]*o[2] + m[1]*o[6] + m[2]*o[10] + m[3]*o[14]
	r[3] = m[0]*o[3] + m[1]*o[7] + m[2]*o[11] + m[3]*o[15]

	r[4] = m[4]*o[0] + m[5]*o[4] + m[6]*o[8] + m[7]*o[12]
	r[5] = m[4]*o[1] + m[5]*o[5] + m[6]*o[9] + m[7]*o[13]
	r[6] = m[4]*o[2] + m[5]*o[6] + m[6]*o[10] + m[7]*o[14]
	r[7] = m[4]*o[3] + m[5]*o[7] + m[6]*o[11] + m[7]*o[15]

	r[8] = m[8]*o[0] + m[9]*o[4] + m[10]*o[8] + m[11]*o[12]
	r[9] = m[8]*o[1] + m[9]*o[5] + m[10]*o[9] + m[11]*o[13]
	r[10] = m[8]*o[2] + m[9]*o[6] + m[10]*o[10] + m[11]*o[14]
	r[11] = m[8]*o[3] + m[9]*o[7] + m[10]*o[11] + m[11]*o[15]

	r[12] = m[12]*o[0] + m[13]*o[4] + m[14]*o[8] + m[15]*o[12]
	r[13] = m[12]*o[1] + m[13]*o[5] + m[14]*o[9] + m[15]*o[13]
	r[14] = m[12]*o[2] + m[13]*o[6] + m[14]*o[10] + m[15]*o[14]
	r[15] = m[12]*o[3] + m[13]*o[7] + m[14]*o[11] + m[15]*o[15]
}

//------------------------------------------------------------------------------

// `LookAt` returns a transform from world space into the specific eye space
// that the projective matrix functions (Perspective, OrthographicFrustum, ...)
// are designed to expect.
//
// See also `Perspective` and `OrthographicFrustum`.
func LookAt(eye, center, up Vec3) Mat4 {
	center.Subtract(eye)
	f := center.Normalized()
	u := up.Normalized()
	s := f.Cross(u).Normalized()
	u = s.Cross(f)

	res := MakeMat4(
		s.X, s.Y, s.Z, -s.Dot(eye),
		u.X, u.Y, u.Z, -u.Dot(eye),
		-f.X, -f.Y, -f.Z, f.Dot(eye),
		0, 0, 0, 1,
	)

	return res
}
