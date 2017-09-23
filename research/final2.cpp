#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <random>

using namespace std;

const uint32_t mod = 11;

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
	
	friend ringe& operator += (ringe &lhs, const ringe &rhs) { return lhs = lhs + rhs; }
	friend ringe& operator -= (ringe &lhs, const ringe &rhs) { return lhs = lhs - rhs; }
	friend ringe& operator *= (ringe &lhs, const ringe &rhs) { return lhs = lhs * rhs; }
	
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
		
	friend constexpr ringe inverse(const ringe &p) {
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
	
	friend constexpr ringe pow(ringe x, uint32_t n){
		ringe y(1);
		for(; n; n /= 2) {
			if (n & 1) y = x * y;
			x = x * x;
		}
		return y;
	}

	friend constexpr ringe sqrt(ringe a) {
		if (mod % 4 == 3)
			return pow(a, (mod + 1) / 4);
		assert(mod % 4 == 1);
		auto inv = inverse(a);
		const uint32_t lmb = (mod - 1) & (1 - mod);
		const uint32_t t = (mod - 1) / lmb;
		auto x = pow(a, (t + 1) / 2);
		ringe c = pow(ringe::first_nr(), t);
		for (uint32_t s = lmb / 4; s; s /= 2){
			auto d = x * x * inv;
			if (d == 1)
				return x;
			if (pow(d, s) == mod - 1)
				x *= c;
			c *= c;
		}
		return x;
	}

	friend constexpr int legendre(ringe a) { return pow(a, (mod - 1) / 2) == 1 ? 1 : -1; }

	friend ostream& operator << (ostream &out, ringe &p) { return out << p.x_; }
private:
	uint32_t x_;
	
	static constexpr ringe first_nr() {
		ringe r(2);
		for (; legendre(r) > 0; ++r.x_);
		return r;
	}
};

template <int n> struct xternion {
public:
	enum { dim = n + 1 };
	xternion() {}
	template <class X> xternion(initializer_list<X> list) {
		int i = 0;
		for (const auto &x : list) {
			assert(i < dim);
			c[i++] = x;
		}
	}
	ringe<mod> operator[] (int i) const { assert(0 <= i && i < dim); return c[i]; }
	ringe<mod>& operator[] (int i) { assert(0 <= i && i < dim); return c[i]; }
	bool operator == (const xternion &q) { return !(*this != q); }
	bool operator != (const xternion &q) {
		for (int i = 0; i < dim; ++i)
			if (c[i] != q.c[i])
				return true;
		return false;
	}
	friend ostream& operator << (ostream &out, const xternion &q) {
		out << '(';
		for (int i = 0; i < dim; ++i) {
			if (i) out << ' ';
			out << q[i];
		}
		return out << ')';
	}
	xternion& operator += (const xternion &rhs) {
		for (int i = 0; i < dim; ++i)
			c[i] = c[i] + rhs[i];
		return *this;
	}
	xternion& operator -= (const xternion &rhs) {
		for (int i = 0; i < dim; ++i)
			c[i] = c[i] - rhs[i];
		return *this;
	}
	friend xternion operator + (xternion lhs, const xternion &rhs) {
		lhs += rhs;
		return lhs;
	}
	friend xternion operator - (xternion lhs, const xternion &rhs) {
		lhs -= rhs;
		return lhs;
	}
	friend xternion operator * (const xternion &lhs, const xternion &rhs) {
		xternion res;
		res[0] = lhs[0] * rhs[0];
		for (int i = 1; i <= n; ++i) {
			for (int j = 1; j <= n; ++j) {
				res[0] += lhs[i] * lu.t[i][j] * rhs[j];
			}
			res[i] += lhs[i] * rhs[0] + lhs[0] * rhs[i];
		}
		return res;
	}
	friend xternion zero_inv(const xternion &lhs) {
		// ringe<mod> a[dim][dim];
		// for (int i = 0; i < dim; ++i)
			// for (int j = 0; j < dim; ++j)
				// a[i ^ j][j] = lhs[i] * lu.t[i & j];
		// for (int j = 0; j < dim; ++j) {
			// for (int i = j; ; ++i) {
				// if (i == dim) {
					// xternion r;
					// for (i = 0; i < j; ++i)
						// r[i] = -a[i][j];
					// r[j] = 1;
					// return r;
				// }
				// if (a[i][j] != 0) {
					// if (i != j)
						// for (int k = j; k < dim; ++k)
							// swap(a[i][k], a[j][k]);
					// break;
				// }
			// }
			// if (a[j][j] != 1) {
				// auto inv = inverse(a[j][j]);
				// for (int i = j; i < dim; ++i)
					// a[j][i] *= inv;
				// assert(a[j][j] == 1);
			// }
			// for (int i = 0; i < dim; ++i)
				// if (i != j) {
					// auto t = a[i][j];
					// for (int k = j; k < dim; ++k)
						// a[i][k] -= a[j][k] * t;
					// assert(a[i][j] == 0);
				// }
		// }
		return xternion();
	}
	static vector<pair<xternion, xternion>> zero_divisors() {
		vector<pair<xternion, xternion>> zds;
		for (int i = 1; i <= n; ++i)
			for (int j = 1; j < i; ++j) {
				auto t = sqrt(lu.t[i][i] * inverse(lu.t[j][j]));
				xternion le, ri;
				le[i] = ri[i] = 1;
				le[j] = t;
				ri[j] = -t;
				zds.emplace_back(le, ri);
			}
		return zds;
	}
private:
	ringe<mod> c[dim];
	
	struct lut_t {
		ringe<mod> t[dim][dim];
		lut_t() {
			int i = 0;
			for (uint32_t x = 2; i < n; ++x)
				if (legendre(ringe<mod>(x)) < 0) {
					++i;
					t[i][i] = x;
					for (int j = 1; j < i; ++j) {
						auto z = sqrt(t[i][i] * t[j][j]);
						t[i][j] = t[j][i] = z;
					}
				}
			for (int i = 1; i <= n; ++i)
				for (int j = 1; j <= n; ++j) {
					assert(t[i][j] * t[j][i] == t[i][i] * t[j][j]);
				}
		}
	};
	
	static lut_t lu;
};

template <int n> typename xternion<n>::lut_t xternion<n>::lu;

const int n = 4;

mt19937_64 rng;
uniform_int_distribution<uint32_t> rnd(0, mod - 1);

xternion<n> xrnd() {
	xternion<n> r;
	uniform_int_distribution<uint32_t> rndi(0, xternion<n>::dim - 1);
	for (int k = 0; k < 2; ++k) {
		r[rndi(rng)] = rnd(rng);
	}
	return r;
}

int main() {
	auto zds = xternion<n>::zero_divisors();
	for (auto zd : zds) {
		assert(zd.first * zd.second == xternion<n>());
	}
	for (int i = 0; i < 1000000; ++i) {
		auto a = xrnd();
		auto b = xrnd();
		auto c = xrnd();
		//cout << a << ' ' << b << ' ' << c << endl;
		assert(a * b == b * a);
		//assert((a * b) * c == a * (b * c));
		assert((a + b) * c == a * c + b * c);
		for (const auto &zd : zds) {
			const auto &l = zd.first;
			const auto &r = zd.second;
			cout << l * (r * a) << endl;
			assert(l * r * a == xternion<n>());
		}
		// auto d = zero_inv(a);
		// if (d != xternion<n>()) {
			// //a = zero_inv(d);
			// cout << a << endl;
			// cout << d << endl;
			// cout << a * d << endl;
			// cout << endl;
		// }
	}
	return 0;
}