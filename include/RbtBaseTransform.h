/***********************************************************************
* The rDock program was developed from 1998 - 2006 by the software team 
* at RiboTargets (subsequently Vernalis (R&D) Ltd).
* In 2006, the software was licensed to the University of York for 
* maintenance and distribution.
* In 2012, Vernalis and the University of York agreed to release the 
* program as Open Source software.
* This version is licensed under GNU-LGPL version 3.0 with support from
* the University of Barcelona.
* http://rdock.sourceforge.net/
***********************************************************************/

//Abstract base class for all classes which manipulate molecular coords

#ifndef _RBTBASETRANSFORM_H_
#define _RBTBASETRANSFORM_H_

#include "RbtConfig.h"
#include "RbtBaseObject.h"

#if 1 // AdvanceSoft(2018)
extern int m_iSdf;
extern int m_iIsland;
extern int counter_iCycle;
extern bool b_islandGA_w_migration;
extern bool initial_b_islandGA_w_migration;
#endif // AdvanceSoft (2018)

class RbtTransformAgg;//forward declaration

class RbtBaseTransform : public RbtBaseObject
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names

  ////////////////////////////////////////
  //Constructors/destructors
  virtual ~RbtBaseTransform();
  
  //Give aggregates access to RbtBaseSF private methods/data
  friend class RbtTransformAgg;

  ////////////////////////////////////////
  //Public methods
  ////////////////
  //Fully qualified name, prefixed by all ancestors
  RbtString GetFullName() const;

  //Main public method - actually apply the transform
  //Not virtual. Base class method checks if transform is enabled,
  //sends SFRequests to reconfigure the scoring function
  //then calls private virtual method Execute() to apply the transform
  void Go();
  
  //Aggregate handling methods
  virtual void Add(RbtBaseTransform*) throw (RbtError);
  virtual void Remove(RbtBaseTransform*) throw (RbtError);
  virtual RbtBool isAgg() const;
  virtual RbtUInt GetNumTransforms() const;
  virtual RbtBaseTransform* GetTransform(RbtUInt) const throw (RbtError);
  void Orphan();//Force removal from the parent aggregate
  RbtBaseTransform* GetParentTransform() const;
	
  //Scoring function request handling
  //Transforms can store up a list of requests to send to the workspace
  //scoring function each time Go() is executed
  void AddSFRequest(RbtRequestPtr);
  void ClearSFRequests();
  void SendSFRequests();

#if 1 // AdvanceSoft (2018)
  static void Set_m_iSdf(int iSdf) { m_iSdf = iSdf; }
  int Get_m_iSdf() { return m_iSdf; }
  static void Set_m_iIsland(int iIsland) { m_iIsland = iIsland; }
  int Get_m_iIsland() { return m_iIsland; }
#endif // AdvanceSoft(2018)
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  RbtBaseTransform(const RbtString& strClass, const RbtString& strName);

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtBaseTransform();
  RbtBaseTransform(const RbtBaseTransform&);//Copy constructor disabled by default
  RbtBaseTransform& operator=(const RbtBaseTransform&);//Copy assignment disabled by default
  
  //PURE VIRTUAL - subclasses should override. Applies the transform
  virtual void Execute() = 0;
  
 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtBaseTransform* m_parent;
  RbtRequestList m_SFRequests;




};

//Useful typedefs
typedef vector<RbtBaseTransform*> RbtBaseTransformList;//Vector of smart pointers
typedef RbtBaseTransformList::iterator RbtBaseTransformListIter;
typedef RbtBaseTransformList::const_iterator RbtBaseTransformListConstIter;

#endif //_RBTBASETRANSFORM_H_
