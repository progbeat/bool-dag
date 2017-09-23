#include <type_traits>

namespace crypto {

template <class byte_t = uint8_t> class datastream
{
public:	
	template <class X>
	typename std::enable_if<std::is_integral<X>::value, datastream&>::type
		operator << (const X &x)
	{
		for (auto i = sizeof x; i--; )
			data_.push_back(static_cast<byte_t>(x >> i * 8));
		return *this;
	}
	template <class X>
	typename std::enable_if<!std::is_integral<X>::value, datastream&>::type
		operator << (const X &x)
	{
		data_.push_back(x);
	}
	const byte_t* data() const { return data_.data(); }
	size_t size() const { return data_.size(); }
private:
	std::vector<byte_t> data_;
};

} // namespace crypto