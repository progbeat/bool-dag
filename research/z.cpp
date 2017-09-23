#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <random>

using namespace std;

#include "ring.h"

const uint32_t mod = 17;

template <uint32_t mod> class ringe {
public:
	template <class X>
	static ringe make(X x) {
		x %= (X)mod;
		if (x < X()) x += mod;
		return ringe((uint32_t)x);
	}
	ringe() : x_(0) {}
	ringe(uint32_t x) {
		assert(x < mod);
		x_ = x;
	}
	ringe(const ringe &rhs) : x_(rhs.x_) {}
	ringe& operator = (const ringe &rhs) {
		x_ = rhs.x_;
		return *this;
	}
	friend ringe operator * (const ringe &p, const ringe &q) {
		return ringe((uint32_t)((uint64_t)p.x_ * q.x_ % mod));
	}
	friend ringe operator + (const ringe &p, const ringe &q) {
		return ringe(p.x_ < mod - q.x_ ? p.x_ + q.x_ : p.x_ - (mod - q.x_));
	}
	friend ringe operator - (const ringe &p, const ringe &q) {
		return ringe(p.x_ < q.x_ ? p.x_ + (mod - q.x_) : p.x_ - q.x_);
	}
	ringe operator - () const { return make(x_ ? mod - x_ : uint32_t()); }
	
	operator uint32_t() const { return x_; }
	
	bool operator == (const ringe &rhs) const { return x_ == rhs.x_; }
	
	bool operator != (const ringe &rhs) const { return x_ != rhs.x_; }
	
	template <class X> bool operator == (X x) const { return x == x_; }
	
	template <class X> bool operator != (X x) const { return x != x_; }
		
	friend ringe inverse(const ringe &p) {
		uint32_t x = p.x_, y = mod;
		uint32_t a = 1, b = 0;
		for (; x; ) {
			uint32_t t = y / x;
			y %= x;
			b -= a * t;
			std::swap(x, y);
			std::swap(a, b);
		}
		return make(y == 1 ? a == mod ? b + mod : b : uint32_t());
	}
	
	static const ringe zero;
	
	friend ostream& operator << (ostream &out, ringe &p) { return out << p.x_; }
private:
	uint32_t x_;
};

const auto a = ringe<mod>::make(-1);
const auto b = ringe<mod>::make(-1);
const ringe<mod> ab = a * b;

struct quaternion {
public:
	quaternion() {}
	quaternion(	const ringe<mod> &c0,
				const ringe<mod> &c1,
				const ringe<mod> &c2,
				const ringe<mod> &c3)
		: c({c0, c1, c2, c3})
	{}
	ringe<mod> operator[] (int i) const {
		assert(0 <= i && i < 4);
		return c[i];
	}
	ringe<mod>& operator[] (int i) {
		assert(0 <= i && i < 4);
		return c[i];
	}
	bool operator != (const quaternion &q) {
		for (int i = 0; i < 4; ++i)
			if (c[i] != q.c[i])
				return true;
		return false;
	}
	bool operator == (const quaternion &q) { return !(*this != q); }
	friend ostream& operator << (ostream &out, const quaternion &q) {
		return out << '(' << q[0] << " + " << q[1] << "i + " << q[2] << "j + " << q[3] << "k)";
	}
	static quaternion zero, one;
private:
	ringe<mod> c[4];
};

quaternion quaternion::zero;
quaternion quaternion::one {1, 0, 0, 0};

quaternion operator * (const quaternion &p, const quaternion &q) {
	quaternion r;
	r[0] = p[0] * q[0] + p[1] * q[1] * a + p[2] * q[2] * b - p[3] * q[3] * ab;
	r[1] = p[0] * q[1] + p[1] * q[0] - p[2] * q[3] * b + p[3] * q[2] * b;
	r[2] = p[0] * q[2] + p[2] * q[0] + p[1] * q[3] * a - p[3] * q[1] * a;
	r[3] = p[0] * q[3] + p[3] * q[0] + p[1] * q[2] - p[2] * q[1];
	return r;
}

quaternion operator + (const quaternion &p, const quaternion &q) {
	return quaternion {
		p[0] + q[0],
		p[1] + q[1],
		p[2] + q[2],
		p[3] + q[3]
	};
}

quaternion operator - (const quaternion &p, const quaternion &q) {
	return quaternion {
		p[0] - q[0],
		p[1] - q[1],
		p[2] - q[2],
		p[3] - q[3]
	};
}

quaternion operator - (const quaternion &p) {
	return quaternion{-p[0], -p[1], -p[2], -p[3]};
}

quaternion l0codiv(const quaternion &q) {
	// 0 = x0 * q[0] + x1 * q[1] * a + x2 * q[2] * b - x3 * q[3] * ab;
	// 0 = x0 * q[1] + x1 * q[0]     - x2 * q[3] * b + x3 * q[2] * b;
	// 0 = x0 * q[2] + x2 * q[0]     + x1 * q[3] * a - x3 * q[1] * a;
	// 0 = x0 * q[3] + x3 * q[0]     + x1 * q[2]     - x2 * q[1];
	ringe<mod> f[4][4] = {
		{q[0], a * q[1],  b * q[2], -ab * q[3]},
		{q[1],     q[0], -b * q[3],   b * q[2]},
		{q[2], a * q[3],      q[0],  -a * q[1]},
		{q[3],     q[2],     -q[1],       q[0]}
	};
	for (int j = 0; j < 4; ++j) {
		int i = j;
		for (; i < 4 && f[i][j] == 0; ++i);
		if (i == 4) {
			quaternion r;
			r[j] = 1;
			for (i = 0; i < j; ++i)
				r[i] = -f[i][j];
			return r;
		}
		if (i != j) {
			for (int k = j; k < 4; ++k)
				swap(f[i][k], f[j][k]);
		}
		if (f[j][j] != 1) {
			auto t = inverse(f[j][j]);
			for (int k = j; k < 4; ++k)
				f[j][k] = f[j][k] * t;
		}
		for (i = 0; i < 4; ++i) if (i != j) {
			auto t = f[i][j];
			for (int k = j; k < 4; ++k)
				f[i][k] = f[i][k] - f[j][k] * t;
		}
	}
	return quaternion();
}

quaternion inverse(const quaternion &q) {
	auto i = inverse(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
	quaternion r {q[0], -q[1], -q[2], -q[3]};
	return quaternion{r[0] * i, r[1] * i, r[2] * i, r[3] * i};
}

mt19937_64 rng;
uniform_int_distribution<uint32_t> rnd(0, mod - 1);

quaternion rndq() {
	quaternion r;
	for (int i = 0; i < 4; ++i)
		r[i] = rnd(rng);
	return r;
}

int main() {
	quaternion r, t;
	for (; ; ) {
		t = rndq();
		auto l = l0codiv(t);
		if (t != quaternion::zero)
			break;		
	}
	for (int k = 0; k < 10000; ++k) {
		r = rndq();
		auto l = l0codiv(r);
		if (l != quaternion::zero) {
			cout << "l = " << l << endl;
			cout << "r = " << r << endl;
			cout << "t = " << t << endl;
			cout << "l * r = " << l * r << endl;
			cout << "l * t * r = " << l * t * r << endl;
			cout << endl;
			// cout << "r = " << r << endl;
			// cout << "* = " << l * r << endl;
			// cout << endl;
		}
	}
	return 0;
}
 
 // -ax + by + cz + dt
 //  ay + bx + ct - dz
 //  az - bt + cx + dy
 //  at + bz - cy + dx
 
// вместо того, чтобы пойти спать я взял и доказал эту херню... так что я больше не мудак)
// Система
// -ax+by+cz+dt=0
 // ay+bx+ct-dz=0
 // az-bt+cx+dy=0
 // at+bz-cy+dx=0
// имеет нетривиальное решение, если матрица
// (-x  y  z  t)
// ( y  x  t -z)
// ( z -t  x  y)
// ( t  z -y  x)
// вырожденная, но детерминант такой матрицы -(x^2+y^2+z^2+t^2)^2
 
 
 // (-a)
 // (b, 0, 0 )
 // (c, 0 )
 // ( 0, 0, 0, +d) * x
 // +
 // (0, +b, 0, 0 )
 // ( a, 0, 0, 0 )
 // ( 0, 0, 0, +d)
 // ( 0, 0, -c, 0) * y
 // +
 // ( 0, 0, c, 0 )
 // ( 0, 0, 0, -d)
 // ( a, 0, 0, 0 )
 // ( 0, b, 0, 0 ) * z
 // +
 // (0, 0, 0, d)
 // (0, 0, c, 0)
 // (0,-b, 0, 0)
 // (a, 0, 0, 0) * t
 