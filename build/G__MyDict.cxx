// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__MyDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "/Users/alepisani/Documents/tans/tracker/include/point.h"
#include "/Users/alepisani/Documents/tans/tracker/include/traccia.h"
#include "/Users/alepisani/Documents/tans/tracker/include/tracklet.h"
#include "/Users/alepisani/Documents/tans/tracker/include/evento.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_point(void *p = nullptr);
   static void *newArray_point(Long_t size, void *p);
   static void delete_point(void *p);
   static void deleteArray_point(void *p);
   static void destruct_point(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::point*)
   {
      ::point *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::point >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("point", ::point::Class_Version(), "include/point.h", 7,
                  typeid(::point), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::point::Dictionary, isa_proxy, 4,
                  sizeof(::point) );
      instance.SetNew(&new_point);
      instance.SetNewArray(&newArray_point);
      instance.SetDelete(&delete_point);
      instance.SetDeleteArray(&deleteArray_point);
      instance.SetDestructor(&destruct_point);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::point*)
   {
      return GenerateInitInstanceLocal(static_cast<::point*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::point*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_traccia(void *p = nullptr);
   static void *newArray_traccia(Long_t size, void *p);
   static void delete_traccia(void *p);
   static void deleteArray_traccia(void *p);
   static void destruct_traccia(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::traccia*)
   {
      ::traccia *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::traccia >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("traccia", ::traccia::Class_Version(), "include/traccia.h", 7,
                  typeid(::traccia), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::traccia::Dictionary, isa_proxy, 4,
                  sizeof(::traccia) );
      instance.SetNew(&new_traccia);
      instance.SetNewArray(&newArray_traccia);
      instance.SetDelete(&delete_traccia);
      instance.SetDeleteArray(&deleteArray_traccia);
      instance.SetDestructor(&destruct_traccia);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::traccia*)
   {
      return GenerateInitInstanceLocal(static_cast<::traccia*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::traccia*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_tracklet(void *p = nullptr);
   static void *newArray_tracklet(Long_t size, void *p);
   static void delete_tracklet(void *p);
   static void deleteArray_tracklet(void *p);
   static void destruct_tracklet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::tracklet*)
   {
      ::tracklet *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::tracklet >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("tracklet", ::tracklet::Class_Version(), "include/tracklet.h", 8,
                  typeid(::tracklet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::tracklet::Dictionary, isa_proxy, 4,
                  sizeof(::tracklet) );
      instance.SetNew(&new_tracklet);
      instance.SetNewArray(&newArray_tracklet);
      instance.SetDelete(&delete_tracklet);
      instance.SetDeleteArray(&deleteArray_tracklet);
      instance.SetDestructor(&destruct_tracklet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::tracklet*)
   {
      return GenerateInitInstanceLocal(static_cast<::tracklet*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::tracklet*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_evento(void *p = nullptr);
   static void *newArray_evento(Long_t size, void *p);
   static void delete_evento(void *p);
   static void deleteArray_evento(void *p);
   static void destruct_evento(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::evento*)
   {
      ::evento *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::evento >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("evento", ::evento::Class_Version(), "include/evento.h", 7,
                  typeid(::evento), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::evento::Dictionary, isa_proxy, 4,
                  sizeof(::evento) );
      instance.SetNew(&new_evento);
      instance.SetNewArray(&newArray_evento);
      instance.SetDelete(&delete_evento);
      instance.SetDeleteArray(&deleteArray_evento);
      instance.SetDestructor(&destruct_evento);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::evento*)
   {
      return GenerateInitInstanceLocal(static_cast<::evento*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::evento*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr point::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *point::Class_Name()
{
   return "point";
}

//______________________________________________________________________________
const char *point::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::point*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int point::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::point*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *point::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::point*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *point::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::point*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr traccia::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *traccia::Class_Name()
{
   return "traccia";
}

//______________________________________________________________________________
const char *traccia::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::traccia*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int traccia::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::traccia*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *traccia::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::traccia*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *traccia::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::traccia*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr tracklet::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *tracklet::Class_Name()
{
   return "tracklet";
}

//______________________________________________________________________________
const char *tracklet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::tracklet*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int tracklet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::tracklet*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *tracklet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::tracklet*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *tracklet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::tracklet*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr evento::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *evento::Class_Name()
{
   return "evento";
}

//______________________________________________________________________________
const char *evento::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::evento*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int evento::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::evento*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *evento::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::evento*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *evento::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::evento*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void point::Streamer(TBuffer &R__b)
{
   // Stream an object of class point.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(point::Class(),this);
   } else {
      R__b.WriteClassBuffer(point::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_point(void *p) {
      return  p ? new(p) ::point : new ::point;
   }
   static void *newArray_point(Long_t nElements, void *p) {
      return p ? new(p) ::point[nElements] : new ::point[nElements];
   }
   // Wrapper around operator delete
   static void delete_point(void *p) {
      delete (static_cast<::point*>(p));
   }
   static void deleteArray_point(void *p) {
      delete [] (static_cast<::point*>(p));
   }
   static void destruct_point(void *p) {
      typedef ::point current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::point

//______________________________________________________________________________
void traccia::Streamer(TBuffer &R__b)
{
   // Stream an object of class traccia.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(traccia::Class(),this);
   } else {
      R__b.WriteClassBuffer(traccia::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_traccia(void *p) {
      return  p ? new(p) ::traccia : new ::traccia;
   }
   static void *newArray_traccia(Long_t nElements, void *p) {
      return p ? new(p) ::traccia[nElements] : new ::traccia[nElements];
   }
   // Wrapper around operator delete
   static void delete_traccia(void *p) {
      delete (static_cast<::traccia*>(p));
   }
   static void deleteArray_traccia(void *p) {
      delete [] (static_cast<::traccia*>(p));
   }
   static void destruct_traccia(void *p) {
      typedef ::traccia current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::traccia

//______________________________________________________________________________
void tracklet::Streamer(TBuffer &R__b)
{
   // Stream an object of class tracklet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(tracklet::Class(),this);
   } else {
      R__b.WriteClassBuffer(tracklet::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_tracklet(void *p) {
      return  p ? new(p) ::tracklet : new ::tracklet;
   }
   static void *newArray_tracklet(Long_t nElements, void *p) {
      return p ? new(p) ::tracklet[nElements] : new ::tracklet[nElements];
   }
   // Wrapper around operator delete
   static void delete_tracklet(void *p) {
      delete (static_cast<::tracklet*>(p));
   }
   static void deleteArray_tracklet(void *p) {
      delete [] (static_cast<::tracklet*>(p));
   }
   static void destruct_tracklet(void *p) {
      typedef ::tracklet current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::tracklet

//______________________________________________________________________________
void evento::Streamer(TBuffer &R__b)
{
   // Stream an object of class evento.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(evento::Class(),this);
   } else {
      R__b.WriteClassBuffer(evento::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_evento(void *p) {
      return  p ? new(p) ::evento : new ::evento;
   }
   static void *newArray_evento(Long_t nElements, void *p) {
      return p ? new(p) ::evento[nElements] : new ::evento[nElements];
   }
   // Wrapper around operator delete
   static void delete_evento(void *p) {
      delete (static_cast<::evento*>(p));
   }
   static void deleteArray_evento(void *p) {
      delete [] (static_cast<::evento*>(p));
   }
   static void destruct_evento(void *p) {
      typedef ::evento current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::evento

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_libMyDict() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_libMyDict) = RecordReadRules_libMyDict();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_libMyDict));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_libMyDict_Impl() {
    static const char* headers[] = {
"/Users/alepisani/Documents/tans/tracker/include/point.h",
"/Users/alepisani/Documents/tans/tracker/include/traccia.h",
"/Users/alepisani/Documents/tans/tracker/include/tracklet.h",
"/Users/alepisani/Documents/tans/tracker/include/evento.h",
nullptr
    };
    static const char* includePaths[] = {
"/opt/homebrew/Cellar/root/6.36.04_1/include/root",
"/Users/alepisani/Documents/tans/tracker/src",
"/Users/alepisani/Documents/tans/tracker",
"/opt/homebrew/Cellar/root/6.36.04_1/include/root",
"/Users/alepisani/Documents/tans/tracker/build/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libMyDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$/Users/alepisani/Documents/tans/tracker/include/point.h")))  point;
class __attribute__((annotate("$clingAutoload$/Users/alepisani/Documents/tans/tracker/include/traccia.h")))  traccia;
class __attribute__((annotate("$clingAutoload$/Users/alepisani/Documents/tans/tracker/include/tracklet.h")))  tracklet;
class __attribute__((annotate("$clingAutoload$/Users/alepisani/Documents/tans/tracker/include/evento.h")))  evento;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libMyDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "/Users/alepisani/Documents/tans/tracker/include/point.h"
#include "/Users/alepisani/Documents/tans/tracker/include/traccia.h"
#include "/Users/alepisani/Documents/tans/tracker/include/tracklet.h"
#include "/Users/alepisani/Documents/tans/tracker/include/evento.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"evento", payloadCode, "@",
"point", payloadCode, "@",
"traccia", payloadCode, "@",
"tracklet", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libMyDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libMyDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libMyDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libMyDict() {
  TriggerDictionaryInitialization_libMyDict_Impl();
}
