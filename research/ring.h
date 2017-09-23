static uint32_t gcd(uint32_t x, uint32_t y) {
	for (; y; ) {
		x %= y;
		if (x == 0)
			return y;
		y %= x;
	}
	return x;
}

template <uint32_t mod> struct zring {

	static uint32_t mul(uint32_t x, uint32_t y) {
		return (uint64_t)x * y % mod;
	}

	static uint32_t sub(uint32_t x, uint32_t y) {
		return x < y ? x + (mod - y) : x - y;
	}

	static uint32_t add(uint32_t x, uint32_t y) {
		return x < mod - y ? x + y : x - (mod - y);
	}

	static uint32_t egcd(uint32_t x, uint32_t y, uint32_t &x1, uint32_t &y1, uint32_t &x0, uint32_t &y0) {
		int64_t ax = 1, bx = 0, ay = 0, by = 1;
		for (; x; ) {
			uint32_t t = y / x;
			y %= x;
			ay -= ax * t;
			by -= bx * t;
			std::swap(x, y);
			std::swap(ax, ay);
			std::swap(bx, by);
		}
		auto nrm = [](int64_t x) -> uint32_t {
			if (x < 0)
				x += mod;
			assert(x < mod);
			return static_cast<uint32_t>(x);
		};
		x0 = nrm(ax);
		y0 = nrm(bx);
		x1 = nrm(ay);
		y1 = nrm(by);
		return y;
	}

	static uint32_t gcdmod(uint32_t x, uint32_t &inv) {
		uint32_t y = mod;
		int64_t a = 1, b = 0;
		for (; x; ) {
			uint32_t t = y / x;
			y %= x;
			b -= a * t;
			std::swap(x, y);
			std::swap(a, b);
		}
		inv = b < 0 ? b + mod : b;
		return y;
	}
	
	static uint32_t inv(uint32_t x) {
		uint32_t r;
		return gcdmod(x, r) == 1 ? r : 0;
	}
	
	static uint32_t negate(uint32_t x) {
		return x ? mod - x : 0;
	}
	
	static void pmulc(uint32_t *r, const uint32_t *p, uint32_t c, size_t n) {
		for (size_t i = 0; i < n; ++i)
			r[i] = mul(p[i], c);
	}

	static void paddmulc(uint32_t *r, const uint32_t *p, uint32_t c, size_t n) {
		for (size_t i = 0; i < n; ++i) {
			r[i] = add(r[i], mul(p[i], c));
		}
	}

	static void psubmulc(uint32_t *r, const uint32_t *p, uint32_t c, size_t n) {
		for (size_t i = 0; i < n; ++i)
			r[i] = sub(r[i], mul(p[i], c));
	}

	static void padd(uint32_t *r, const uint32_t *p, const uint32_t *q, size_t n) {
		for (size_t i = 0; i < n; ++i)
			r[i] = add(p[i], q[i]);
	}

	static void psub(uint32_t *r, const uint32_t *p, const uint32_t *q, size_t n) {
		for (size_t i = 0; i < n; ++i)
			r[i] = sub(p[i], q[i]);
	}
	
	static void pmul22(uint32_t *x, uint32_t *y, uint32_t x1, uint32_t y1, uint32_t x0, uint32_t y0, size_t n) {
		for (size_t i = 0; i < n; ++i) {
			const uint32_t xi = x[i], yi = y[i];
			x[i] = add(mul(xi, x1), mul(yi, y1));
			y[i] = add(mul(xi, x0), mul(yi, y0));
		}
	}
	
	static int geliminate(uint32_t *a, int n, int m, int k = 0) {
		const int stride = m;
		for (; m > k; --m, --n, a += stride + 1) {
			for (int i = 0; ; ++i) {
				if (i == n)
					return stride - m;
				uint32_t *b = a + i * stride;
				if (b[0]) {
					if (i) {
						for (int j = 0; j < m; ++j)
							std::swap(a[j], b[j]);
					}
					break;
				}
			}
			if (uint32_t t = inv(a[0])) {
				pmulc(a + 1, a + 1, t, m - 1);
				a[0] = 1;
			}
			for (int i = 1; i < n; ++i) {
				uint32_t *b = a + i * stride;
				if (b[0] == 0)
					continue;
				if (a[0] == 1) {
					uint32_t b0 = b[0];
					b[0] = 0;
					psubmulc(b + 1, a + 1, b0, m - 1);
				} else {
					uint32_t cj1, ci1, cj0, ci0;
					uint32_t c, d, e, f, g = egcd(a[0], b[0], c, d, e, f);
					assert(add(mul(a[0], c), mul(b[0], d)) == g);
					assert(add(mul(a[0], e), mul(b[0], f)) == 0);
					a[0] = g;
					b[0] = 0;
					pmul22(a + 1, b + 1, c, d, e, f, m - 1);
				}
			}
		}
		return stride - m;
	}

	static void null_vector(uint32_t *r, uint32_t *a, int d, int m) {
		assert(d < m);
		r[d] = 1;
		uint32_t *c = &a[d];
		for (int i = d; i --> 0; ) {
			uint32_t *b = &a[i];
			if (b[i * m] == 1) {
				r[i] = negate(c[i * m]);
				for (int j = 0; j < i; ++j) {
					c[j * m] = add(c[j * m], mul(b[j * m], r[i]));
				}
			} else {
				uint32_t bf, cf, t[2];
				egcd(b[i * m], c[i * m], t[0], t[1], bf, cf);
				for (int j = 0; j < i; ++j) {
					c[j * m] = add(mul(b[j * m], bf), mul(c[j * m], cf));
				}
				r[i] = bf;
				pmulc(&r[i + 1], &r[i + 1], cf, d - i);
			}
		}
	}
	
private:
	zring() {}
};

template <uint32_t Mod> class poly : public vector<uint32_t> {
public:
	typedef zring<Mod> ring;

	poly() {}

	template <class X> poly(const std::initializer_list<X> &list) {
		reserve(list.size());
		for (auto x : list) {
			if (x < 0)
				x = ring::negate((-x) % Mod);
			else
				x %= Mod;
			push_back(static_cast<uint32_t>(x < 0 ? x + Mod : x));
		}
		norm();
	}

	template <class X> poly(const X &other)
		: vector<uint32_t>(other.begin(), other.end())
	{
		norm();
	}

	friend poly& operator += (poly &p, const poly &q) {
		if (p.size() < q.size())
			p.resize(q.size());
		ring::padd(&p[0], &p[0], &q[0], q.size());
		return p.norm();
	}

	friend poly& operator -= (poly &p, const poly &q) {
		if (p.size() < q.size())
			p.resize(q.size());
		ring::psub(&p[0], &p[0], &q[0], q.size());
		return p.norm();
	}

	friend poly operator + (const poly &p, const poly &q) {
		poly r(p);
		r += q;
		return r;
	}

	friend poly operator - (const poly &p, const poly &q) {
		poly r(p);
		r -= p;
		return r;
	}

	friend poly operator * (const poly &p, const poly &q) {
		if (p.empty() || q.empty())
			return poly();
		poly r(p.size() + q.size() - 1);
		for (size_t i = 0; i < p.size(); ++i)
			if (p[i]) {
				ring::paddmulc(&r[i], &q[0], p[i], q.size());
			}
		r.norm();
		return r;
	}

	friend poly& operator *= (poly &p, const poly &q) { return p = p * q; }

	friend poly divmod(poly &p, const poly &q) {
		assert(!q.empty() && q.back());
		if (p.size() < q.size())
			return poly();
		poly d(p.size() - q.size() + 1);
		const size_t qn = q.size() - 1;
		if (q.back() == 1) {
			for (size_t i = p.size(); i -->= q.size(); )
				if (p[i]) {
					d[i - qn]  = p[i];
					ring::psubmulc(&p[i - qn], &q[0], p[i]);
				}
		} else {
			uint32_t inv = ring::inv(q.back());
			for (size_t i = p.size(); i -->= q.size(); )
				if (p[i]) {
					uint32_t c = mul(p[i], inv);
					d[i - qn] = c;
					ring::psubmulc(&p[i - qn], &q[0], c);
				}
		}
		p.norm();
		d.norm();
		return d;
	}

	friend poly& operator %= (poly &p, const poly &q) {
		assert(!q.empty() && q.back());
		if (p.size() < q.size())
			return p;
		if (q.back() == 1) {
			for (size_t i = p.size(); i -->= q.size(); )
				if (p[i]) {
					ring::psubmulc(&p[i - (q.size() - 1)], &q[0], p[i], q.size() - 1);
					p[i] = 0;
				}
		} else {
			uint32_t inv = ring::inv(q.back());
			assert(inv != 0);
			for (size_t i = p.size(); i -->= q.size(); )
				if (p[i]) {
					uint32_t c = ring::mul(p[i], inv);
					ring::psubmulc(&p[i + 1 - q.size()], &q[0], c, q.size());
					assert(p[i] == 0);
				}
		}
		return p.norm();
	}

	friend poly operator / (const poly &p, const poly &q) {
		poly t(p);
		return divmod(t, q);
	}

	friend poly& operator /= (poly &p, const poly &q) { return p = p / q; }

	friend poly operator % (const poly &p, const poly &q) {
		poly t(p);
		t %= q;
		return t;
	}

	friend ostream& operator << (ostream &out, const poly &p) {
		out << '(';
		for (size_t i = 0; i < p.size(); ++i) {
			if (i)
				out << ' ';
			out << p[i];
		}
		return out << ')';
	}

	poly null(const poly &mod) {
		if (empty() || mod.empty())
			return poly();
		if (mod.back() != 1) {
			uint32_t inv = ring::inv(mod.back());
			if (inv == 0)
				return poly();
			poly nmod = mod;
			ring::pmulc(&nmod[0], &mod[0], inv, mod.size());
			return null(nmod);
		}
		const int n = mod.degree();
		std::vector<uint32_t> a(n * n);
		{
			poly p = *this;
			p.resize(n);
			for (int i = 0; ; ) {
				uint32_t *po = &a[i];
				for (int j = 0; j < n; ++j)
					po[j * n] = p[j];
				if (++i == n)
					break;
				uint32_t mf = p[n - 1];
				std::move_backward(&p[0], &p[n - 1], &p[n]);
				p[0] = 0;
				ring::psubmulc(&p[0], &mod[0], mf, n);
			}
		}
		int d = ring::geliminate(&a[0], n, n);
		if (d == n)
			return poly();
		poly r(static_cast<size_t>(d + 1));
		ring::null_vector(&r[0], &a[0], d, n);
		return r;
	}

	int degree() const {
		size_t n = size();
		return static_cast<int>(n ? n - 1 : 0);
	}
	
	template <class Rng> poly random_zero_divisor(Rng &rng) const;
	
	uint32_t operator()(uint32_t x) const {
		uint32_t y = 0;
		for (auto i = rbegin(); i != rend(); ++i) {
			y = ring::mul(x, y);
			y = ring::add(*i, y);
		}
		return y;
	}
	
private:
	poly(size_t deg) : vector<uint32_t>(deg) {}

	poly& norm() { 
		for (; !empty() && !back(); pop_back());
		return *this;
	}
};

template <uint32_t Mod> template <class Rng>
poly<Mod> poly<Mod>::random_zero_divisor(Rng &rng) const {
	if (empty())
		return poly<Mod>();
	assert(back() == 1);
	const int n = degree(), n1 = n + 1;
	std::vector<uint32_t> a(n * n), b(n * n1), f(n), g(n);
	f[n / 2] = 1;
	for (int i = 1; i < n / 2; ++i) {
		f[i] = rng();
		assert(f[i] < Mod);
	}
	std::copy(f.begin(), f.end(), a.begin());
	for (int i = 1; i < n; ++i) {
		uint32_t *t = &a[i * n];
		t[0] = 0;
		std::copy(t - n, t - 1, t + 1);
		ring::psubmulc(t, &front(), t[-1], n);
	}
	// for (int i = 0; i < n; ++i) {
		// for (int j = 0; j < n; ++j) {
			// cout << a[i * n + j] << ' ';
		// }
		// cout << endl;
	// }
	// cout << endl;
	for (int i = 0; i < n; ++i) {
		b[i * n1] = rng();
		assert(b[i * n1] < Mod);
	}
	for (int k = 0; k < n; ++k) {
		std::fill(g.begin(), g.end(), 0);
		uint32_t *c = &b[k];
		for (int i = 0; i < n; ++i) {
			uint32_t ci = c[i * n1];
			uint32_t *ai = &a[i * n];
			for (int j = 0; j < n; ++j) {
				g[j] = ring::add(g[j], ring::mul(ai[j], ci));
			}
		}
		++c;
		for (int i = 0; i < n; ++i)
			c[i * n1] = g[i];
	}
	if (ring::geliminate(&b[0], n, n1) != n)
		return poly();
	poly r((size_t)n1);
	ring::null_vector(&r[0], &b[0], n, n1);
	//cout << r << endl;
	std::vector<uint32_t> roots;
	for (uint32_t lambda = 0; lambda < Mod; ++lambda) {
		if (r(lambda) == 0) {
			roots.push_back(lambda);
		}
	}
	cout << "cp: " << r << endl;
	if (roots.size()) {
		cout << "roots:";
		for (auto x : roots)
			cout << ' ' << x;
		cout << endl;
		f[0] = ring::negate(roots[0]);
		return f;
	}
	return poly();
}