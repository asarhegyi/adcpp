#ifndef SPLIT_T
#define SPLIT_T

#include <cstddef>

struct skip_empty
{
	enum empties_t { empties_ok, no_empties };
};


template <typename Container>
Container& split(
	Container&                            result,
	const typename Container::value_type& str,
	const typename Container::value_type& delimiters,
	skip_empty::empties_t                 empties = skip_empty::empties_ok )
{
	result.clear();
	size_t current;
	size_t next = -1;
	do
	{
		if (empties == skip_empty::no_empties)
		{
			next = str.find_first_not_of( delimiters, next + 1 );
			if (next == Container::value_type::npos) break;
			next -= 1;
		}
		current = next + 1;
		next = str.find_first_of( delimiters, current );
		result.push_back( str.substr( current, next - current ) );
	}
	while (next != Container::value_type::npos);

	return result;
}

#endif
