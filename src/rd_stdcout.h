#pragma once
#include <iostream>

template < typename C, typename T = std::char_traits<C> >
struct basic_teebuf : public std::basic_streambuf<C, T>
{
	typedef std::basic_streambuf<C, T> streambuf_type;
	typedef typename T::int_type int_type;

	basic_teebuf(streambuf_type* buff_a, streambuf_type* buff_b)
		: first(buff_a), second(buff_b) {}

protected:
	virtual int_type overflow(int_type c)
	{
		const int_type eof = T::eof();
		if (T::eq_int_type(c, eof)) return T::not_eof(c);
		else
		{
			const C ch = T::to_char_type(c);
			if (T::eq_int_type(first->sputc(ch), eof) ||
				T::eq_int_type(second->sputc(ch), eof))
				return eof;
			else return c;
		}
	}

	virtual int sync()
	{
		return !first->pubsync() && !second->pubsync() ? 0 : -1;
	}

private:
	streambuf_type* first;
	streambuf_type* second;
};

template < typename C, typename T = std::char_traits<C> >
struct basic_teestream : public std::basic_ostream<C, T>
{
	typedef std::basic_ostream<C, T> stream_type;
	typedef basic_teebuf<C, T> streambuff_type;

	basic_teestream(stream_type& first, stream_type& second)
		: stream_type(&stmbuf), stmbuf(first.rdbuf(), second.rdbuf()) {}

	basic_teestream(streambuff_type* first, streambuff_type* second)
		: stream_type(&stmbuf), stmbuf(first, second) {}

	~basic_teestream() { stmbuf.pubsync(); }

private: streambuff_type stmbuf;
};

typedef basic_teebuf<char> teebuf;

typedef basic_teestream<char> teestream;



/// <summary>
/// Source: https://cplusplus.com/forum/general/97966/#msg532935
/// </summary>
class rdCoutStream
{
public:
	std::streambuf* stdoutbuf = nullptr;
	std::filebuf* fbuf = nullptr;
	teebuf* tbuf = nullptr;
public:
	rdCoutStream(std::string file) {
		fbuf = new std::filebuf;
		fbuf->open(file, std::ios::out);
		stdoutbuf = std::cout.rdbuf();
		tbuf=new teebuf(fbuf, stdoutbuf);
		std::cout.rdbuf(tbuf);
	}

	~rdCoutStream() {
		std::cout.rdbuf(stdoutbuf);
		delete tbuf;
		delete fbuf;
	}
};

