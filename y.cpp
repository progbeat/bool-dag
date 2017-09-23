#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <random>

using namespace std;

const uint32_t mod = 0x00c7818b;
//const uint32_t mod = 105647;

uint32_t mulx(uint32_t x, uint32_t z = mod)
{
	return (x << 1) ^ ((int32_t(x) >> 31) & z);
}

uint32_t mul(uint32_t x, uint32_t y, uint32_t z = mod)
{
	uint32_t r = 0;
	for (int i = 0; i < 32; ++i)
	{
		r ^= x & -(y & 1);
		y >>= 1;
		x = mulx(x);
	}
	return r;
}

pair<uint32_t, uint32_t> div(uint32_t z, uint32_t x)
{
	uint32_t q = 0;
	int i = __builtin_clz(x) + 1, j = 31;
	q |= uint32_t(1) << i;
	z ^= x << i;
	for (; i --> 0; --j)
		if (z & uint32_t(1) << j) {
			q |= uint32_t(1) << i;
			z ^= x << i;
		}
	return make_pair(q, z);
}

uint32_t dmod(uint32_t z, uint32_t x)
{
	int i = __builtin_clz(x), j = 31;
	for (; i >= 0; --i, --j)
		if (z & uint32_t(1) << j)
			z ^= x << i;
	assert(z < x);
	return z;
}

bool is_prime(uint32_t z)
{
	if ((z & 1) == 0)
		return false;
	uint32_t up = z;
	for (uint32_t x = 3; x < up; x += 2) {
		if (dmod(z, x) == 0)
			return false;
	}
	return true;
}

int main()
{
	vector<uint32_t> primes;
	vector<pair<uint32_t, uint32_t>> vs;
	for (uint32_t x = 3; x < 0x20000; x += 2) {
		auto d = div(mod, x);
		if (d.second == 0 && x < d.first) {
			if (is_prime(x))
				primes.push_back(x);
			if (is_prime(d.first))
				primes.push_back(d.first);				
			vs.emplace_back(x, d.first);
		}
	}
	uint32_t y = 1;
	for (uint32_t x : primes) {
		cout << x << ' ';
		y = mul(x, y);
		cout << y << ' ';
		y = mul(x, y);
		cout << y << ' ';
		y = mul(x, y);
		cout << y << ' ';
		y = mul(x, y);
		cout << y << ' ';
		y = mul(x, y);
		cout << y << ' ';
		cout << endl;
	}
	cout << primes.size() << endl;
	return 0;
	const int n = (int)vs.size();
	mt19937_64 rng;
	double A = 0;
	int B = 0;
	for (; ; ) {
		++B;
		for (int i = 0; i < n; ++i) {
			int j = uniform_int_distribution<int>(0, i * 2 + 1)(rng);
			if (j > i) {
				swap(vs[i].first, vs[i].second);
				j -= i + 1;
			}
			assert(i >= j);
			swap(vs[i], vs[j]);
		}
		uint32_t x = 1;
		for (int i = 0; i < n; ++i) {
			++A;
			x = mul(x, vs[i].first);
			if (x == 0)
				break;
		}
		cout << A / B << endl;
	}
	return 0;
}