#include <algorithm>
#include <iostream>
#include <cassert>
#include <vector>
#include <random>
#include <ctime>
#include <set>

using namespace std;

#include "ring.h"

//const uint32_t mod = 2 * 3 * 5 * 7 * 11 * 13;
const uint32_t mod = 210 * 47 * 53 * 257;

//mt19937_64 rng(time(0));
mt19937_64 rng;

poly<mod> rnd_poly(size_t n, bool coprime = true) {
	std::vector<uint32_t> p(n);
	std::uniform_int_distribution<uint32_t> rnd(0, mod - 1);
	for (; ; ) {
		uint32_t g = 0;
		for (size_t i = 0; i < n; ++i) {
			do p[i] = rnd(rng);
			while (coprime && zring<mod>::inv(p[i]) == 0);
			g = gcd(g, p[i]);
		}
		if (n < 2 || g == 1)
			return p;
	}
}

pair<double, int> polyqual(const poly<mod> &pm) {
	const int tries = 1000000;
	const int dim = 30;
	int chunk = 0, ra = 0, rb = 0, tot = 0;
	int best_chunk = 0;
	std::uniform_int_distribution<uint32_t> rnd(0, (1 << dim) - 1);
	poly<mod> pr[dim];
	for (int i = 0; i < dim; ++i)
		pr[i] = {1};
	for (int t = 0; t < tries; ++t) {
		poly<mod> p = rnd_poly(pm.degree(), false);
		auto q = p.null(pm);
		if (!q.empty()) {
			++tot;
			if (!((p * q) % pm).empty()) {
				cout << "pm = " << pm << endl;
				cout << "p = " << p << endl;
				cout << "q = " << q << endl;
				cout << ((p * q) % pm) << endl;
				exit(1);
			}
			bool good = true;
			poly<mod> npr[dim];
			uint32_t msk = rnd(rng);
			for (int j = 0; j < dim; ++j) {
				npr[j] = (pr[j] * ((msk >> j) & 1 ? p : q)) % pm;
				if (npr[j].empty()) {
					good = false;
					break;
				}
			}
			if (good) {
				for (int j = 0; j < dim; ++j) {
					pr[j].swap(npr[j]);
				}
				++chunk;
			} else {
				best_chunk = max(best_chunk, chunk);
				ra += chunk;
				++rb;
				chunk = 0;
				for (int j = 0; j < dim; ++j)
					pr[j] = {1};
			}
		}
	}
	return make_pair(best_chunk, tot);
}

void doit(const poly<mod> &pm) {
	const int tries = 10000000;
	const int degree = pm.degree();
	set<poly<mod>> all;
	vector<poly<mod>> ps, qs;
	for (int t = 0; t < tries; ++t) {
		poly<mod> p = rnd_poly(degree, false);
		auto q = p.null(pm);
		if (!q.empty() && !all.count(q)) {
			p = q.null(pm);
			q = p.null(pm);
			if (p != q && !all.count(p) && !all.count(q)) {
				assert((p * q % pm).empty());
				all.insert(p);
				all.insert(q);
				ps.push_back(p);
				qs.push_back(q);
				cout << ps.size() << ' ';
			}
		}
	}
	cout << endl;
	const int dim = 30, n = (int)ps.size();
	std::uniform_int_distribution<uint32_t> rnd(0, (1 << dim) - 1);
	
	struct ff { poly<mod> f[dim]; };
	vector<ff> f[n + 1];
	f[0].emplace_back();
	for (int i = 0; i < dim; ++i)
		f[0][0].f[i] = {1};
	for (int i = 0; i < n; ++i) {
		uint32_t msk = rnd(rng);
		for (int j = i; j >= 0; --j) {
			for (ff &t : f[j]) {
				ff nt;
				bool ok = true;
				for (int k = 0; k < dim; ++k) {
					nt.f[k] = t.f[k] * (msk & 1u << k ? ps[i] : qs[i]) % pm;
					if (nt.f[k].empty()) {
						ok = false;
						break;
					}
				}
				if (ok) {
					f[j + 1].push_back(nt);
				}
			}
			if (f[j + 1].size())
				cout << j + 1 << ": " << f[j + 1].size() << endl;
		}
	}
}

void find_good_mod() {
	pair<double, int> bqual;
	auto rnd = []() -> uint32_t { return rng() % mod; };
	for (int i = 1; i < mod; ++i) {
	#if 1
		poly<mod> pm = rnd_poly(5);
		pm.push_back(1);
	#else
		poly<mod> pm = {i, 0, 0, 0, 1};
	#endif
		// cout << "pmod = " << pm << endl;
		// auto p = pm.random_zero_divisor(rnd);
		// if (p.empty())
			// continue;
		// auto q = p.null(pm);
		// if (!q.empty()) {
			// cout << "ok!" << endl;
			// cout << p << endl;
			// cout << q << endl;
			// p = q.null(pm);
			// cout << p << endl;
			// q = p.null(pm);
			// cout << q << endl;
			// cout << endl;
			// assert((p * q % pm).empty());
		// }
		// else {
			// cout << "fuck" << endl;
			// cout << pm << endl;
			// cout << p << endl;
			// cout << q << endl;
			// return;
		// }
		auto qual = polyqual(pm);
		if (qual > bqual) {
			bqual = qual;
		}
		cout << pm << ' ' << qual.second << ' ' << qual.first << ' ' << bqual.first << endl;
		//return;
	}
}

int main() {
	cout << "mod = " << mod << endl;
	cout.precision(4);
	cout.setf(ios::fixed);
	poly<mod> pm = {1};
	pm *= poly<mod>{11, 1};
	pm *= poly<mod>{13, 1};
	pm *= poly<mod>{17, 1};
	pm *= poly<mod>{19, 1};
	pm *= poly<mod>{23, 1};
	pm *= poly<mod>{29, 1};
	pm *= poly<mod>{31, 1};
	doit(pm);
	//find_good_mod();
	return 0;
}