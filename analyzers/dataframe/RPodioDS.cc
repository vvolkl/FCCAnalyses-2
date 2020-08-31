#include <ROOT/RDF/Utils.hxx>
#include <ROOT/TSeq.hxx>
#include "RPodioDS.h"
#include <ROOT/RMakeUnique.hxx>
#include <TError.h>

namespace ROOT {

namespace RDF {

std::vector<void *> RPodioDS::GetColumnReadersImpl(std::string_view, const std::type_info &ti)
{
   // We know we have only one column and that it's holding ULong64_t's
   if (ti != typeid(ULong64_t)) {
      throw std::runtime_error("The type specified for the column \"col0\" is not ULong64_t.");
   }
   std::vector<void *> ret;
   for (auto i : ROOT::TSeqU(fNSlots)) {
      fCounterAddr[i] = &fCounter[i];
      ret.emplace_back((void *)(&fCounterAddr[i]));
   }
   return ret;
}

RPodioDS::RPodioDS(ULong64_t size, bool skipEvenEntries) : fSize(size), fSkipEvenEntries(skipEvenEntries)
{
}

RPodioDS::~RPodioDS()
{
}

const std::vector<std::string> &RPodioDS::GetColumnNames() const
{
   return fColNames;
}

bool RPodioDS::HasColumn(std::string_view colName) const
{
   return colName == fColNames[0];
}

std::string RPodioDS::GetTypeName(std::string_view) const
{
   return "ULong64_t";
}

std::vector<std::pair<ULong64_t, ULong64_t>> RPodioDS::GetEntryRanges()
{
   auto ranges(std::move(fEntryRanges)); // empty fEntryRanges
   return ranges;
}

bool RPodioDS::SetEntry(unsigned int slot, ULong64_t entry)
{
   if (fSkipEvenEntries && 0 == entry % 2) {
      return false;
   }
   fCounter[slot] = entry;
   return true;
}

void RPodioDS::SetNSlots(unsigned int nSlots)
{
   R__ASSERT(0U == fNSlots && "Setting the number of slots even if the number of slots is different from zero.");

   fNSlots = nSlots;
   fCounter.resize(fNSlots);
   fCounterAddr.resize(fNSlots);
}

void RPodioDS::Initialise()
{
   const auto chunkSize = fSize / fNSlots;
   auto start = 0UL;
   auto end = 0UL;
   for (auto i : ROOT::TSeqUL(fNSlots)) {
      start = end;
      end += chunkSize;
      fEntryRanges.emplace_back(start, end);
      (void)i;
   }
   // TODO: redistribute reminder to all slots
   fEntryRanges.back().second += fSize % fNSlots;
}

std::string RPodioDS::GetLabel()
{
   return "TrivialDS";
}

RInterface<RDFDetail::RLoopManager, RPodioDS> MakePodioDataFrame(ULong64_t size, bool skipEvenEntries)
{
   auto lm = std::make_unique<RDFDetail::RLoopManager>(std::make_unique<RPodioDS>(size, skipEvenEntries),
                                                       RDFInternal::ColumnNames_t{});
   return RInterface<RDFDetail::RLoopManager, RPodioDS>(std::move(lm));
}

} // ns RDF

} // ns ROOT

int main() {

  return 0;
  }
