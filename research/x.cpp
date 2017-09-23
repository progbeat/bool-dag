#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cassert>

using namespace std;

uint64_t mul(uint32_t x, uint32_t y)
{
	uint64_t r = 0;
	for (; y; ) {
		uint64_t i = y & -y;
		y &= y - 1;
		r ^= x * i;
	}
	return r;
}

int main()
{
	vector<short> cnt(0x80000000);
	for (uint32_t l = 0x80000000, x = 3; x < 0x20000; x += 2) {
		cout << x << ' ';
		for (uint64_t z = uint64_t(x) * l; z >= 0x200000000; ) {
			z >>= 1;
			l >>= 1;
		}
		for (uint32_t y = max(l | 1, x); y != 1; y += 2) {
			uint64_t z = mul(x, y);
			if (z < 0x200000000) {
				assert(z >= 0x100000000);
				cnt[(z - 0x100000000) / 2]++;
			} else {
				break;
			}
		}
	}
	int best = 0;
	for (uint32_t i = 0; i < 0x80000000; ++i) {
		if (cnt[i] > best) {
			best = cnt[i];
			cout << i * 2 + 1 << ' ' << best << endl;
		}
	}
	return 0;
}