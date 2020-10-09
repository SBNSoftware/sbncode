////////////////////////////////////////////////////////////////////////////
// \file  AssociationUtil.h
// \brief Utility object to perform functions of association 
//        Ported from the NOvA Utilities package
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////////

#ifndef ASSOCIATIONUTIL_H
#define ASSOCIATIONUTIL_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/EDFilter.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Provenance/ProductID.h" 
#include "canvas/Persistency/Common/FindOne.h" 
#include "canvas/Persistency/Common/FindOneP.h" 

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace util {

  /// \name Create and query associations
  /// see https://cdcvs.fnal.gov/redmine/projects/art/wiki/Inter-Product_References
  /// for information about using art::Assns
  //@{

  /// \brief Create a 1 to 1 association between a new product and one already in the event
  ///
  /// \param a The collection about to be added to the event
  /// \param b The product already in the event
  /// \param indx Which element of \a a to associate. By default the last one.
  /// \param instance instance label for product, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDProducer const& prod,
                                                    art::Event            & evt,
                                                    std::vector<T>          & a,
                                                    art::Ptr<U>               b,
                                                    art::Assns<T,U>      & assn,
                                                    size_t        indx=UINT_MAX,
                                                    std::string const& instance=std::string());

  /// \brief Create a 1 to 1 association between each of a series of new
  ///        products and one already in the event
  ///
  /// \param a The collection about to be added to the event
  /// \param b The product already in the event
  /// \param begin_indx Which element of \a a to associate first.
  /// \param end_indx One more than the index of the last element to associate
  /// \param instance instance label for product, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDProducer const& prod,
                                                    art::Event             &evt,
                                                    std::vector<T>           &a,
                                                    art::Ptr<U>               b,
                                                    art::Assns<T,U>       &assn,
                                                    size_t           begin_indx,
                                                    size_t             end_indx,
                                                    std::string const& instance=std::string());

  /// \brief Create a 1 to 1 association between two products already in the event
  ///
  /// \param a A product already in the event
  /// \param b Another product already in the event
  template<class T, class U> static bool CreateAssn(art::EDProducer const& prod,
                                                    art::Event             &evt,
                                                    art::Ptr<T>              &a,
                                                    art::Ptr<U>               b,
                                                    art::Assns<T,U>       &assn);

  /// \brief Create a 1 to many association between a new product and a PtrVector already in the event
  ///
  /// \param a The collection about to be added to the event
  /// \param b The products already in the event
  /// \param indx Which element of \a a to associate. By default the last one.
  /// \param instance instance label for product, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDProducer const& prod,
                                                    art::Event             &evt,
                                                    std::vector<T>           &a,
                                                    art::PtrVector<U>         b,
                                                    art::Assns<T,U>       &assn,
                                                    size_t        indx=UINT_MAX,
                                                    std::string const& instance=std::string());

  /// \brief Create a 1 to many association between products already in the event
  ///
  /// \param a A product already in the event
  /// \param b A vector of products already in the event (the many)
  template<class T, class U> static bool CreateAssn(art::EDProducer const&  prod,
                                                    art::Event              &evt,
                                                    art::Ptr<T>               &a,
                                                    std::vector< art::Ptr<U> > b,
                                                    art::Assns<T,U>        &assn);

  /// \brief Create a 1 to many association between a new product and a vector of Ptrs already in the event
  ///
  /// \param a The collection about to be added to the event
  /// \param b The products already in the event
  /// \param indx Which element of \a a to associate. By default the last one.
  /// \param instance instance label for product, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDProducer const&  prod,
                                                    art::Event              &evt,
                                                    std::vector<T>            &a,
                                                    std::vector< art::Ptr<U> > b,
                                                    art::Assns<T,U>        &assn,
                                                    size_t         indx=UINT_MAX,
                                                    std::string const& instance=std::string());

  /// \brief Create a 1 to many association between new products
  ///
  /// \param a A collection about to be added to the event
  /// \param b Another collection about to be added to the event
  /// \param startU The first element of \a b to associate
  /// \param endU The last element of \a b to associate +1 (like STL begin() and end())
  /// \param indx Which element of \a a to associate. By default the last one.
  /// \param instancea instance label for product a, defaulted to be an empty string
  /// \param instanceb instance label for product b, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDProducer const&  prod,
                                                    art::Event              &evt,
                                                    std::vector<T>            &a,
                                                    std::vector<U>            &b,
                                                    art::Assns<T,U>        &assn,
                                                    size_t                startU,
                                                    size_t                  endU,
                                                    size_t         indx=UINT_MAX,
                                                    std::string const& instancea=std::string(),
                                                    std::string const& instanceb=std::string());

  /// \brief Create a 1 to 1 between new products
  ///
  /// \param a A collection about to be added to the event
  /// \param b Another collection about to be added to the event
  /// \param indxb Which element of \a b to associate.  By default the last one.
  /// \param indxa Which element of \a a to associate. By default the last one.
  /// \param instancea instance label for product a, defaulted to be an empty string
  /// \param instanceb instance label for product b, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDProducer const& prod,
                                                    art::Event             &evt,
                                                    art::Assns<T,U>       &assn,
                                                    std::vector<T>           &a,
                                                    std::vector<U>           &b,
                                                    size_t       indxb=UINT_MAX,
                                                    size_t       indxa=UINT_MAX,
                                                    std::string const& instancea=std::string(),
                                                    std::string const& instanceb=std::string());
  
  /// \brief Return all objects of type U that are not associated to objects of
  /// type T.
  ///
  /// Label is the module label that would have produced
  /// the associations and likely the objects of type T
  /// this method assumes there is a one to many relationship between T and U
  /// for example if you want to get all rb::CellHits
  /// that are not associated to rb::Clusters
  /// std::vector<const rb::CellHit*> hits = FindUNotAssociatedToU<rb::Cluster>(art::Handle<rb::CellHit>, ...);
  template<class T, class U> static std::vector<const U*> FindUNotAssociatedToT(art::Handle<U>     b, 
                                                                                art::Event  const& evt,
                                                                                std::string const& label);

  /// \brief Return all objects of type U that are not associated to objects of
  /// type T.
  ///
  /// Label is the module label that would have produced the associations and
  /// likely the objects of type T this method assumes there is a one to many
  /// relationship between T and U for example if you want to get all rb::CellHits
  /// that are not associated to rb::Clusters std::vector<art::Ptr<rb::CellHit> > hits =
  /// FindUNotAssociatedToTP<rb::Cluster>(art::Handle<rb::CellHit>, ...);
  template<class T, class U> static std::vector< art::Ptr<U> > FindUNotAssociatedToTP(art::Handle<U>         b, 
    art::Event  const&   evt,
    std::string const& label);

  //@}


  /// \name Create and query associations
  /// see https://cdcvs.fnal.gov/redmine/projects/art/wiki/Inter-Product_References
  /// for information about using art::Assns
  //@{

  /// \brief Create a 1 to 1 association between a new product and one already in the event
  ///
  /// \param a The collection about to be added to the event
  /// \param b The product already in the event
  /// \param indx Which element of \a a to associate. By default the last one.
  /// \param instance instance label for product, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDFilter const&    prod,
                                                    art::Event              &evt,
                                                    std::vector<T>            &a,
                                                    art::Ptr<U>                b,
                                                    art::Assns<T,U>        &assn,
                                                    size_t         indx=UINT_MAX,
                                                    std::string const& instance=std::string());

  /// \brief Create a 1 to 1 association between each of a series of new
  ///        products and one already in the event
  ///
  /// \param a The collection about to be added to the event
  /// \param b The product already in the event
  /// \param begin_indx Which element of \a a to associate first.
  /// \param end_indx One more than the index of the last element to associate
  /// \param instance instance label for product, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDFilter const&  prod,
                                                    art::Event            &evt,
                                                    std::vector<T>          &a,
                                                    art::Ptr<U>              b,
                                                    art::Assns<T,U>      &assn,
                                                    size_t          begin_indx,
                                                    size_t            end_indx,
                                                    std::string const& instance=std::string());

  /// \brief Create a 1 to 1 association between two products already in the event
  ///
  /// \param a A product already in the event
  /// \param b Another product already in the event
  template<class T, class U> static bool CreateAssn(art::EDFilter const&  prod,
                                                    art::Event            &evt,
                                                    art::Ptr<T>             &a,
                                                    art::Ptr<U>              b,
                                                    art::Assns<T,U>      &assn);

  /// \brief Create a 1 to many association between a new product and a PtrVector already in the event
  ///
  /// \param a The collection about to be added to the event
  /// \param b The products already in the event
  /// \param indx Which element of \a a to associate. By default the last one.
  /// \param instance instance label for product, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDFilter const&   prod,
                                                    art::Event             &evt,
                                                    std::vector<T>           &a,
                                                    art::PtrVector<U>         b,
                                                    art::Assns<T,U>       &assn,
                                                    size_t        indx=UINT_MAX,
                                                    std::string const& instance=std::string());

  /// \brief Create a 1 to many association between products already in the event
  ///
  /// \param a A product already in the event
  /// \param b A vector of products already in the event (the many)
  template<class T, class U> static bool CreateAssn(art::EDFilter const&    prod,
                                                    art::Event              &evt,
                                                    art::Ptr<T>               &a,
                                                    std::vector< art::Ptr<U> > b,
                                                    art::Assns<T,U>        &assn);

  /// \brief Create a 1 to many association between a new product and a vector of Ptrs already in the event
  ///
  /// \param a The collection about to be added to the event
  /// \param b The products already in the event
  /// \param indx Which element of \a a to associate. By default the last one.
  /// \param instance instance label for product, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDFilter const&    prod,
                                                    art::Event              &evt,
                                                    std::vector<T>            &a,
                                                    std::vector< art::Ptr<U> > b,
                                                    art::Assns<T,U>        &assn,
                                                    size_t         indx=UINT_MAX,
                                                    std::string const& instance=std::string());

  /// \brief Create a 1 to many association between new products
  ///
  /// \param a A collection about to be added to the event
  /// \param b Another collection about to be added to the event
  /// \param startU The first element of \a b to associate
  /// \param endU The last element of \a b to associate +1 (like STL begin() and end())
  /// \param indx Which element of \a a to associate. By default the last one.
  /// \param instancea instance label for product a, defaulted to be an empty string
  /// \param instanceb instance label for product b, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDFilter const&   prod,
                                                    art::Event             &evt,
                                                    std::vector<T>           &a,
                                                    std::vector<U>           &b,
                                                    art::Assns<T,U>       &assn,
                                                    size_t               startU,
                                                    size_t                 endU,
                                                    size_t        indx=UINT_MAX,
                                                    std::string const& instancea=std::string(),
                                                    std::string const& instanceb=std::string());

  /// \brief Create a 1 to 1 between new products
  ///
  /// \param a A collection about to be added to the event
  /// \param b Another collection about to be added to the event
  /// \param indxb Which element of \a b to associate.  By default the last one.
  /// \param indxa Which element of \a a to associate. By default the last one.
  /// \param instancea instance label for product a, defaulted to be an empty string
  /// \param instanceb instance label for product b, defaulted to be an empty string
  template<class T, class U> static bool CreateAssn(art::EDFilter const& prod,
                                                    art::Event           &evt,
                                                    art::Assns<T,U>     &assn,
                                                    std::vector<T>         &a,
                                                    std::vector<U>         &b,
                                                    size_t     indxb=UINT_MAX,
                                                    size_t     indxa=UINT_MAX,
                                                    std::string const& instancea=std::string(),
                                                    std::string const& instanceb=std::string());

}// end namespace

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDProducer const& prod,
                                                        art::Event           & evt,
                                                        std::vector<T>       & a,
                                                        art::Ptr<U>            b,
                                                        art::Assns<T,U>      & assn,
                                                        size_t                 indx,
                                                        std::string     const& instance)
{
  bool ret = true;
  
  if(indx == UINT_MAX) indx = a.size()-1;
  
  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instance);
    art::Ptr<T> aptr(aid, indx, evt.productGetter(aid));
    assn.addSingle(aptr, b);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDProducer const& prod,
                                                        art::Event            &evt,
                                                        std::vector<T>        &a,
                                                        art::Ptr<U>            b,
                                                        art::Assns<T,U>       &assn,
                                                        size_t                 begin_indx,
                                                        size_t                 end_indx,
                                                        std::string    const&  instance)
{
  if(end_indx == UINT_MAX) end_indx = a.size();

  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instance);
    auto getter = evt.productGetter(aid);
    for(size_t i = begin_indx; i < end_indx; ++i){
      art::Ptr<T> aptr(aid, i, getter);
      assn.addSingle(aptr, b);
    }
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    return false;
  }

  return true;
}


//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDProducer const& prod,
                                                        art::Event            &evt,
                                                        art::Ptr<T>           &a,
                                                        art::Ptr<U>            b,
                                                        art::Assns<T,U>       &assn)
{
  bool ret = true;
  
  try{
    assn.addSingle(a, b);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDProducer const& prod,
                                                        art::Event            &evt,
                                                        std::vector<T>        &a,
                                                        art::PtrVector<U>      b,
                                                        art::Assns<T,U>       &assn,
                                                        size_t                 indx,
                                                        std::string    const&  instance)
{
  
  bool ret = true;
  
  if(indx == UINT_MAX) indx = a.size()-1;
  
  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instance);
    art::Ptr<T> aptr(aid, indx, evt.productGetter(aid));
    for(size_t i = 0; i < b.size(); ++i) assn.addSingle(aptr, b[i]);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDProducer const&    prod,
                                                        art::Event               &evt, 
                                                        std::vector<T>           &a,
                                                        std::vector<art::Ptr<U> > b,
                                                        art::Assns<T,U>          &assn,
                                                        size_t                    indx,
                                                        std::string     const&    instance)
{
  bool ret = true;

  if(indx == UINT_MAX) indx = a.size()-1;

  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instance);
    art::Ptr<T> aptr(aid, indx, evt.productGetter(aid));

    for(size_t i = 0; i < b.size(); ++i) assn.addSingle(aptr, b[i]);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
                                      << "art:Assns, exception thrown: "
                                      << e;
    ret = false;
  }

  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDProducer const&    prod,
                                                        art::Event               &evt, 
                                                        art::Ptr<T>              &a,
                                                        std::vector<art::Ptr<U> > b,
                                                        art::Assns<T,U>          &assn)
{
  bool ret = true;

  try{
    for(size_t i = 0; i < b.size(); ++i) assn.addSingle(a, b[i]);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
                                      << "art:Assns, exception thrown: "
                                      << e;
    ret = false;
  }

  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDProducer const& prod,
                                                        art::Event            &evt,
                                                        std::vector<T>        &a,
                                                        std::vector<U>        &b,
                                                        art::Assns<T,U>       &assn,
                                                        size_t                 startU,
                                                        size_t                 endU,
                                                        size_t                 indx,
                                                        std::string     const& instancea,
                                                        std::string     const& instanceb)
{
  
  bool ret = true;
  
  if(indx == UINT_MAX) indx = a.size()-1;
  
  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instancea);
    art::Ptr<T> aptr(aid, indx, evt.productGetter(aid));
    for(size_t i = startU; i < endU; ++i){
      art::ProductID bid = evt.getProductID< std::vector<U> >(instanceb);
      art::Ptr<U> bptr(bid, i, evt.productGetter(bid));
      assn.addSingle(aptr, bptr);
    }
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDProducer const& prod,
                                                        art::Event            &evt,
                                                        art::Assns<T,U>       &assn,
                                                        std::vector<T>        &a,
                                                        std::vector<U>        &b,
                                                        size_t                 indxb,
                                                        size_t                 indxa,
                                                        std::string     const& instancea,
                                                        std::string     const& instanceb)
{
  
  bool ret = true;
  
  if(indxa == UINT_MAX) indxa = a.size()-1;
  if(indxb == UINT_MAX) indxb = b.size()-1;
  
  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instancea);
    art::Ptr<T> aptr(aid, indxa, evt.productGetter(aid));
    art::ProductID bid = evt.getProductID< std::vector<U> >(instanceb);
    art::Ptr<U> bptr(bid, indxb, evt.productGetter(bid));
    assn.addSingle(aptr, bptr);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline std::vector<const U*> util::FindUNotAssociatedToT(art::Handle<U>     b,
                                                                                    art::Event  const& evt,
                                                                                    std::string const& label)
{

  // Do a FindOne for type T for each object of type U
  // If the FindOne returns an invalid maybe ref, add the pointer
  // of object type U to the return vector

  std::vector<const U*> notAssociated;

  art::FindOne<T> fa(b, evt, label);

  for(size_t u = 0; u < b->size(); ++u){
    cet::maybe_ref<T const> t(fa.at(u));
    if( !t.isValid() ){
      art::Ptr<U> ptr(b, u);
      notAssociated.push_back(ptr.get());
    }
  }
  
  return notAssociated;
}

//----------------------------------------------------------------------
template<class T, class U> inline std::vector< art::Ptr<U> > util::FindUNotAssociatedToTP(art::Handle<U>     b,
                                                                                          art::Event  const& evt,
                                                                                          std::string const& label)
{

  // Do a FindOneP for type T for each object of type U
  // If the FindOne returns an invalid maybe ref, add the pointer
  // of object type U to the return vector

  std::vector< art::Ptr<U> > notAssociated;

  art::FindOneP<T> fa(b, evt, label);

  for(size_t u = 0; u < b->size(); ++u){
    cet::maybe_ref<T const> t(fa.at(u));
    if( !t.isValid() ){
      art::Ptr<U> ptr(b, u);
      notAssociated.push_back(ptr);
    }
  }
  
  return notAssociated;
}




//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDFilter const& prod,
                                                        art::Event            &evt,
                                                        std::vector<T>        &a,
                                                        art::Ptr<U>            b,
                                                        art::Assns<T,U>       &assn,
                                                        size_t                 indx,
                                                        std::string     const& instance)
{
  bool ret = true;
  
  if(indx == UINT_MAX) indx = a.size()-1;
  
  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instance);
    art::Ptr<T> aptr(aid, indx, evt.productGetter(aid));
    assn.addSingle(aptr, b);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDFilter const& prod,
                                                        art::Event            &evt,
                                                        std::vector<T>        &a,
                                                        art::Ptr<U>            b,
                                                        art::Assns<T,U>       &assn,
                                                        size_t                 begin_indx,
                                                        size_t                 end_indx,
                                                        std::string    const&  instance)
{
  if(end_indx == UINT_MAX) end_indx = a.size();

  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instance);
    auto getter = evt.productGetter(aid);
    for(size_t i = begin_indx; i < end_indx; ++i){
      art::Ptr<T> aptr(aid, i, getter);
      assn.addSingle(aptr, b);
    }
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    return false;
  }

  return true;
}


//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDFilter const& prod,
                                                        art::Event            &evt,
                                                        art::Ptr<T>           &a,
                                                        art::Ptr<U>            b,
                                                        art::Assns<T,U>       &assn)
{
  bool ret = true;
  
  try{
    assn.addSingle(a, b);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDFilter const& prod,
                                                        art::Event            &evt,
                                                        std::vector<T>        &a,
                                                        art::PtrVector<U>      b,
                                                        art::Assns<T,U>       &assn,
                                                        size_t                 indx,
                                                        std::string    const&  instance)
{
  
  bool ret = true;
  
  if(indx == UINT_MAX) indx = a.size()-1;
  
  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instance);
    art::Ptr<T> aptr(aid, indx, evt.productGetter(aid));
    for(size_t i = 0; i < b.size(); ++i) assn.addSingle(aptr, b[i]);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDFilter const&    prod,
                                                        art::Event               &evt, 
                                                        std::vector<T>           &a,
                                                        std::vector<art::Ptr<U> > b,
                                                        art::Assns<T,U>          &assn,
                                                        size_t                    indx,
                                                        std::string     const&    instance)
{
  bool ret = true;

  if(indx == UINT_MAX) indx = a.size()-1;

  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instance);
    art::Ptr<T> aptr(aid, indx, evt.productGetter(aid));

    for(size_t i = 0; i < b.size(); ++i) assn.addSingle(aptr, b[i]);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
                                      << "art:Assns, exception thrown: "
                                      << e;
    ret = false;
  }

  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDFilter const&    prod,
                                                        art::Event               &evt, 
                                                        art::Ptr<T>              &a,
                                                        std::vector<art::Ptr<U> > b,
                                                        art::Assns<T,U>          &assn)
{
  bool ret = true;

  try{
    for(size_t i = 0; i < b.size(); ++i) assn.addSingle(a, b[i]);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
                                      << "art:Assns, exception thrown: "
                                      << e;
    ret = false;
  }

  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDFilter const& prod,
                                                        art::Event            &evt,
                                                        std::vector<T>        &a,
                                                        std::vector<U>        &b,
                                                        art::Assns<T,U>       &assn,
                                                        size_t                 startU,
                                                        size_t                 endU,
                                                        size_t                 indx,
                                                        std::string     const& instancea,
                                                        std::string     const& instanceb)
{
  
  bool ret = true;
  
  if(indx == UINT_MAX) indx = a.size()-1;
  
  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instancea);
    art::Ptr<T> aptr(aid, indx, evt.productGetter(aid));
    for(size_t i = startU; i < endU; ++i){
      art::ProductID bid = evt.getProductID< std::vector<U> >(instanceb);
      art::Ptr<U> bptr(bid, i, evt.productGetter(bid));
      assn.addSingle(aptr, bptr);
    }
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}

//----------------------------------------------------------------------
template<class T, class U> inline bool util::CreateAssn(art::EDFilter const& prod,
                                                        art::Event           &evt,
                                                        art::Assns<T,U>     &assn,
                                                        std::vector<T>         &a,
                                                        std::vector<U>         &b,
                                                        size_t              indxb,
                                                        size_t              indxa,
                                                        std::string const& instancea,
                                                        std::string const& instanceb)
{
  
  bool ret = true;
  
  if(indxa == UINT_MAX) indxa = a.size()-1;
  if(indxb == UINT_MAX) indxb = b.size()-1;
  
  try{
    art::ProductID aid = evt.getProductID< std::vector<T> >(instancea);
    art::Ptr<T> aptr(aid, indxa, evt.productGetter(aid));
    art::ProductID bid = evt.getProductID< std::vector<U> >(instanceb);
    art::Ptr<U> bptr(bid, indxb, evt.productGetter(bid));
    assn.addSingle(aptr, bptr);
  }
  catch(cet::exception &e){
    mf::LogWarning("AssociationUtil") << "unable to create requested "
				      << "art:Assns, exception thrown: "
				      << e;
    ret = false;
  }
  
  return ret;
}



#endif  //ASSOCIATIONUTIL_H
