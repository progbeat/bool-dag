#include "crypto/sha256.hxx"
#include "crypto/datastream.hxx"

#include <iostream>
#include <cassert>
#include <cstring>
#include <vector>
#include <string>
#include <ctime>

using namespace std;
using namespace crypto;

vector<uint32_t> sha256(const char *str) {
    vector<uint32_t> h(8);
    sha256(h.data(), str, strlen(str));
    return h;
}

vector<uint32_t> sha256(const std::string &str) {
    vector<uint32_t> h(8);
    sha256(h.data(), str.data(), str.size());
    return h;
}

string hex(const vector<uint32_t> &data) {
    string res(data.size() * 8, ' ');
    for (size_t i = 0; i < data.size(); ++i) {
        uint32_t x = data[i];
        for (int j = 8; j --> 0; ) {
            int l = x % 16;
            x /= 16;
            res[i * 8 + j] = (l < 10 ? '0' : 'a' - 10) + l;
        }
    }
    return res;
}

void performance_test() {
    vector<char> msg(0x40000000, 'x');
    vector<uint32_t> h(8);
    clock_t start = clock();
    sha256(&h[0], msg.data(), msg.size());
    clock_t et = clock() - start;
    cout << "Speed = " << CLOCKS_PER_SEC * (1.0 / (1 << 20)) * msg.size() / et << " MiB/s" << endl;
}

void basic_tests()
{
    assert(
        hex(sha256("")) ==
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert(
        hex(sha256("abc")) ==
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    assert(
        hex(sha256("message digest"))
        == "f7846f55cf23e14eebeab5b4e1550cad5b509e3348fbc4efa3a1413d393cb650"
    );
    assert(
        hex(sha256("The quick brown fox jumps over the lazy dog"))
        == "d7a8fbb307d7809469ca9abcb0082e4f8d5651e46d3cdb762d02d0bf37c9e592"
    );
    assert(
        hex(sha256("The quick brown fox jumps over the lazy dog."))
        == "ef537f25c895bfa782526529a9b63d97aa631564d5d789c2b765448c8635fb6c"
    );
    assert(
        hex(sha256("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnop"))
        == "aa353e009edbaebfc6e494c8d847696896cb8b398e0173a4b5c1b636292d87c7"
    );
    assert(
        hex(sha256("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"))
        == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
    );
    assert(
        hex(sha256("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq1"))
        == "cbb143ed5e1ae1ea21653c91cde5c1be208e326ffea9013f98bcea239f214b5b"
    );
    assert(
        hex(sha256("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"))
        == "db4bfcbd4da0cd85a60c3c37d3fbd8805c77f15fc6b1fdfe614ee0a7c8fdb4c0"
    );
    assert(
        hex(sha256("For this sample, this 63-byte string will be used as input data"))
        == "f08a78cbbaee082b052ae0708f32fa1e50c5c421aa772ba5dbb406a2ea6be342"
    );
    assert(
        hex(sha256("This is exactly 64 bytes long, not counting the terminating byte"))
        == "ab64eff7e88e2e46165e29f2bce41826bd4c7b3552f6b382a9e7d3af47c245f8"
    );
    assert(
        hex(sha256("12345678901234567890123456789012345678901234567890123456789012345678901234567890"))
        == "f371bc4a311f2b009eef952dd83ca80e2b60026c8e935592d0f9c308453c813e"
    );
    assert(
        hex(sha256(string(1000000, 'a')))
        == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"
    );
}

void double_hash_tests()
{
	auto h = sha256("hello");
	assert(hex(h) == "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824");
	datastream<uint8_t> ds;
	for (auto x : h) ds << x;
	sha256(&h[0], ds);
	assert(hex(h) == "9595c9df90075148eb06860365df33584b75bff782a510c6cd4883a419833d50");
}

int main()
{
    basic_tests();
	double_hash_tests();
    performance_test();
}