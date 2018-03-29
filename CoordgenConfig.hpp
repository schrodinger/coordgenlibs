#pragma once


#ifdef WIN32
#define EXPORT_COORDGEN __declspec(dllimport)
#else
#define EXPORT_COORDGEN
#endif

