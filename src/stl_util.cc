template<typename InputIterator, typename OutputIterator,typename Predicate>
OutputIterator copy_if(InputIterator begin,
                       InputIterator end,
                       OutputIterator destBegin,
                       Predicate p)
{
  while (begin!=end){
    if(p(*begin))*destBegin++=*begin;
    ++begin;
  }
  return destBegin;
}
