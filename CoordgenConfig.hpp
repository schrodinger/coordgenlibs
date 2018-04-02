#pragma once

#ifdef IN_COORDGEN
#ifdef WIN32
#define EXPORT_COORDGEN __declspec(dllexport)
#else
#define EXPORT_COORDGEN __attribute__((visibility("default")))
#endif

#else

#ifdef WIN32
#define EXPORT_COORDGEN __declspec(dllimport)
#else
#define EXPORT_COORDGEN
#endif

#endif