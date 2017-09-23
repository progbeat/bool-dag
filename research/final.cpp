#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <random>
#include <ctime>

using namespace std;

const uint32_t mod = 1000000007;

mt19937_64 rng(time(0));
uniform_int_distribution<uint32_t> rnd(0, mod - 1);

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
	enum { dim = 1 << n };
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
		for (int i = 0; i < dim; ++i)
			for (int j = 0; j < dim; ++j) {
				auto z = lhs[i] * rhs[i];
				res[i | j] += z * lu.t[i & j].a;
				res[i ^ j] += z * lu.t[i & j].b;
			}
		return res;
	}
	friend xternion zero_inv(const xternion &lhs) {
		ringe<mod> a[dim][dim];
		// for (int i = 0; i < dim; ++i)
			// for (int j = 0; j < dim; ++j)
				// a[i ^ j][j] = lhs[i] * lu.t[i & j];
		for (int j = 0; j < dim; ++j) {
			for (int i = j; ; ++i) {
				if (i == dim) {
					xternion r;
					for (i = 0; i < j; ++i)
						r[i] = -a[i][j];
					r[j] = 1;
					return r;
				}
				if (a[i][j] != 0) {
					if (i != j)
						for (int k = j; k < dim; ++k)
							swap(a[i][k], a[j][k]);
					break;
				}
			}
			if (a[j][j] != 1) {
				auto inv = inverse(a[j][j]);
				for (int i = j; i < dim; ++i)
					a[j][i] *= inv;
				assert(a[j][j] == 1);
			}
			for (int i = 0; i < dim; ++i)
				if (i != j) {
					auto t = a[i][j];
					for (int k = j; k < dim; ++k)
						a[i][k] -= a[j][k] * t;
					assert(a[i][j] == 0);
				}
		}
		return xternion();
	}
	static vector<pair<xternion, xternion>> zero_divisors() {
		vector<pair<xternion, xternion>> zds;
		for (int i = 1; i < dim; ++i) {
			if (!__builtin_parity(i) || (i & i - 1))
				continue;
			//assert(legendre(lu.t[i]) < 0);
			// for (int j = 1; j < i; ++j)
				// if ((i & j) == 0 && __builtin_parity(j) && (j & j - 1) == 0) {
					// auto t = sqrt(lu.t[i] * inverse(lu.t[j]));
					// xternion le, ri;
					// le[i] = ri[i] = 1;
					// le[j] = -t;
					// ri[j] = t;
					// zds.emplace_back(le, ri);
				// }
		}
		return zds;
	}
private:
	ringe<mod> c[dim];
	
	struct lut_t {
		struct { ringe<mod> a, b; } t[dim];
		
		lut_t() {
			t[0] = {0, 1};
			for (int i = 0; i < n; ++i) {
				ringe<mod> a, b, d;
				do {
					a = rnd(rng);
					b = rnd(rng);
					d = a * a + ringe<mod>(4) * b;
				} while (a == 0 || b == 0 || legendre(d) > 0);
				t[1 << i] = {a, b};
			}
			auto inv2 = inverse(ringe<mod>(2));
			for (int j = 1; j < dim; ++j) {
				int i = j & -j;
				if (i == j)
					continue;
				ringe<mod> a1 = t[i].a * inv2;
				ringe<mod> b1 = t[i].b;
				ringe<mod> a2 = t[i ^ j].a * inv2;
				ringe<mod> b2 = t[i ^ j].b;
				ringe<mod> d1 = a1 * a1 + b1;
				ringe<mod> d2 = a2 * a2 + b2;
				assert(legendre(d1) < 0);
				assert(legendre(d2) < 0);
				ringe<mod> f = d1 * d2;
				assert(legendre(f) > 0);
				f = sqrt(f);
				ringe<mod> a1a2 = a1 * a2;
				ringe<mod> a3 = a1a2 + f;
				ringe<mod> d3 = a1 * a1 * d1 + a2 * a2 * d2 + a1a2 * f * ringe<mod>(2);
				assert(legendre(d3) < 0);
				ringe<mod> b3 = d3 - a3 * a3;
				t[j] = {a3 * ringe<mod>(2), b3};
			}
		}
	};
	
	static lut_t lu;
};

template <int n> typename xternion<n>::lut_t xternion<n>::lu;

const int n = 3;

xternion<n> xrnd() {
	xternion<n> r;
	uniform_int_distribution<uint32_t> rndi(0, xternion<n>::dim - 1);
	for (int k = 0; k < 10; ++k) {
		r[rndi(rng)] = rnd(rng);
	}
	return r;
}

// 10 * 9 / 2 = 45
// 11 * 10  55
// 12 * 11 / 2 = 66

int main() {
	// auto zds = xternion<n>::zero_divisors();
	// cout << "#zds = " << zds.size() << endl;
	// for (auto zd : zds) {
		// assert(zd.first * zd.second == xternion<n>());
	// }
	// const int m = zds.size();
	// for (int i = 0; i < m; ++i)
		// for (int j = 0; j < m; ++j) {
			// if (i == j)
				// continue;
			// assert(zds[i].first * zds[j].second != xternion<n>());
		// }
	// const int M = (1 << m);
	// uniform_int_distribution<uint32_t> rnd_msk(0, M - 1);
	// for (int k = 0; ; ++k) {
		// xternion<n> r;
		// r[0] = 1;
		// uint32_t q = rnd_msk(rng);
		// for (int i = 0; i < m; ++i)
			// r = r * zds[i].first;
		// //r = r * r * r * r;
		// if (r == xternion<n>()) {
			// vector<xternion<n>> v;
			// for (int i = 0; i < m; ++i) {
				// cout << zds[i].first << endl;
				// v.push_back(zds[i].first);
			// }
			// cout << endl;
			// for (int i = 0; i < m; ++i) {
				// cout << zds[i].second << endl;
			// }
			// cout << endl;
			// for (int j = 0; j < v.size(); ++j)
				// for (int i = 0; i < j; ++i) {
					// cout << v[i] * v[j] << endl;
				// }
			// cout << endl;
			// cout << r << endl;
		// }
		// assert(r != xternion<n>());
		// // if (k % 1024 == 0) {
			// // cout << k << ' ';
		// // }
	// }
	// return 0;
	for (int i = 0; i < 1000000; ++i) {
		auto a = xrnd();
		auto b = xrnd();
		auto c = xrnd();
		cout << a << endl;
		cout << b << endl;
		cout << c << endl;
		assert(a * b == b * a);
		assert((a * b) * c == a * (b * c));
		assert((a + b) * c == a * c + b * c);
		auto d = zero_inv(a);
		if (d != xternion<n>()) {
			//a = zero_inv(d);
			cout << a << endl;
			cout << d << endl;
			cout << a * d << endl;
			cout << endl;
		}
	}
	return 0;
}