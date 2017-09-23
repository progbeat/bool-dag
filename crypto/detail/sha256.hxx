#include <algorithm>

namespace crypto {
namespace detail {

template <class uint32_t, class uint8_t> class chunks_stream
{
public:
	chunks_stream(const uint8_t *data, size_t n)
		: data_(data)
		, off_(0)
		, n_(n)
	{}
	bool next(uint32_t chunk[])
	{
		if (off_ + 64 <= n_) {
			for (int i = 0; i < 16; ++i)
				chunk[i] =	data_[i * 4 + 0] << 24 |
							data_[i * 4 + 1] << 16 |
							data_[i * 4 + 2] <<  8 |
							data_[i * 4 + 3];
		} else if (off_ <= n_ + 8) {
			std::fill(chunk, chunk + 16, uint32_t());
			if (off_ <= n_) {
				size_t nl = n_ - off_;
				for (size_t i = 0; i < nl; ++i)
					chunk[i / 4] |= data_[i] << (~i & 3) * 8;
				chunk[nl / 4] |= 0x80 << (~nl & 3) * 8;
			}
			if (off_ + 56 > n_) {
				uint64_t len = uint64_t(n_) * 8;
				chunk[14] = len >> 32;
				chunk[15] = len & 0xffffffff;
			}
		}
		else
			return false;
		off_ += 64;
		data_ += 64;
		return true;
	}
private:
	const uint8_t *data_;
	size_t off_;
	size_t n_;
};

template <class X> inline X rots(X x) { return X(); }

template<int i, int... Args, class X> inline X rots(X x)
{
    return (x >> i) ^ (x << (32 - i)) ^ rots<Args...>(x);
}

static uint32_t sha256_k[] = {
	0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
	0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
	0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
	0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
	0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
	0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
	0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
	0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

template <int i> struct sha256_helper
{
	template <class X>
	static void round(X &a, X &b, X &c, X &d, X &e, X &f, X &g, X &h, const X *w)
	{
		X t = h + rots<6, 11, 25>(e) + ((e & f) ^ (~e & g)) + w[i] + sha256_k[i];
		sha256_helper<i + 1>::round(
			h = t + rots<2, 13, 22>(a) + ((a & b) ^ (a & c) ^ (b & c)),
			a, b, c, d += t, e, f, g, w
		);
	}
};

template <> struct sha256_helper<64>
{
	template <class X>
	static void round(X &a, X &b, X &c, X &d, X &e, X &f, X &g, X &h, const X *w) {}
};

} // namespace detail
} // namespace crypto