#ifndef FEDDUTILS_hpp
#define FEDDUTILS_hpp

#include "feddlib/core/FEDDCore.hpp"


#ifndef FEDD_TIMER_START
#define FEDD_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("FEDD") + std::string(S))));
#endif

#ifndef FEDD_TIMER_STOP
#define FEDD_TIMER_STOP(A) A.reset();
#endif


namespace FEDD{
    
template<class ForwardIt, class GO>
ForwardIt uniqueWithCombines(ForwardIt first, ForwardIt last, std::vector<std::vector<GO> >& combines)
{
    
    if (first == last)
        return last;
    
    ForwardIt firstForDistance = first;
    ForwardIt result = first;
    combines[ distance( firstForDistance, result ) ].push_back( (GO) distance( firstForDistance, first ) );
    while (++first != last) {
        if (!(*result == *first) && ++result != first) {
            *result = std::move(*first);
            // also add the element which is the final unique element (the first in the sorted list)
            combines[ distance( firstForDistance, result ) ].push_back( (GO) distance( firstForDistance, first ) );
        }
        else{
            combines[ distance( firstForDistance, result ) ].push_back( (GO) distance( firstForDistance, first ) );
        }
    }
    return ++result;
}

template <typename T>
std::vector<T> sort_from_ref(
                        std::vector<T> const& in,
                        std::vector<int> const& reference
                        ) {
    std::vector<T> ret(in.size());
    
    int const size = in.size();
    for (int i = 0; i < size; ++i)
        ret[i] = in[reference[i]];
    
    return ret;
}

template <typename T>
std::vector<T> sort_from_ref(
                             std::vector<T> const& in,
                             std::vector<long long> const& reference
                             ) {
    std::vector<T> ret(in.size());
    
    int const size = in.size();
    for (long long i = 0; i < size; ++i)
        ret[i] = in[reference[i]];
    
    return ret;
}

    
template <typename T>
void sort2byFirst( std::vector<std::vector<T> >& in, std::vector<T>& in2 )
{

    std::vector<int> index(in.size(), 0);
    for (size_t i = 0 ; i != index.size() ; i++)
        index[i] = i;
    
    std::sort(index.begin(), index.end(),
              [&](const int& a, const int& b) {
                  return  in[a] < in[b];
              }
              );
    in = sort_from_ref( in, index );
    in2 = sort_from_ref( in2, index );
        
}
    
template <typename T, class GO>
void make_unique( std::vector<std::vector<T> >& in, vec2D_GO_Type& combinedElements, std::vector<GO>& globaIDs )
{
    {
        std::vector<int> index(in.size(), 0);
        for (size_t i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        std::sort(index.begin(), index.end(),
             [&](const int& a, const int& b) {
                 return  in[a] < in[b];
             }
             );
        in = sort_from_ref( in, index );
        globaIDs = sort_from_ref( globaIDs, index );
    }
    {
        std::vector<int> index(in.size(), 0);
        for (size_t i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        combinedElements.resize( in.size() );
        
        auto it = uniqueWithCombines( in.begin(), in.end(), combinedElements );
        
        in.resize( distance( in.begin(), it ) );
        combinedElements.resize( in.size() );
    }
}
    
template <typename T>
void make_unique( std::vector<std::vector<T> >& in )
{
    {
        std::vector<int> index(in.size(), 0);
        for (size_t i = 0 ; i != index.size() ; i++)
            index[i] = i;

        std::sort(index.begin(), index.end(),
             [&](const int& a, const int& b) {
                 return  in[a] < in[b];
             }
             );
        in = sort_from_ref( in, index );
    }
    {

        auto it = std::unique( in.begin(), in.end() );

        in.resize( distance( in.begin(), it ) );
    }
}

template <typename T>
void make_unique( std::vector<std::vector<T> >& in, vec2D_GO_Type& combinedElements )
{
    {
        std::vector<int> index(in.size(), 0);
        for (size_t i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        std::sort(index.begin(), index.end(),
             [&](const int& a, const int& b) {
                 return  in[a] < in[b];
             }
             );
        in = sort_from_ref( in, index );
        
    }
    {
        std::vector<int> index(in.size(), 0);
        for (size_t i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        combinedElements.resize( in.size() );
        
        auto it = uniqueWithCombines( in.begin(), in.end(), combinedElements );
        
        in.resize( distance( in.begin(), it ) );
        
        combinedElements.resize( in.size() );
    }
}
    
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());
    
    std::vector<T> result;
    result.reserve(a.size());
    
    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<T>());
    return result;
}
    
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());
    
    std::vector<T> result;
    result.reserve(a.size());
    
    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::minus<T>());
    return result;
}
    
template <typename T>
void make_unique( std::vector<T>& in )
{
    std::sort( in.begin(), in.end() );
    auto it = unique( in.begin(), in.end() );
    in.erase( it, in.end() );
}


}
#endif