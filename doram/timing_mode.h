#pragma once

#ifndef GIGADORAM_ENABLE_VENDOR_TIMING
#define GIGADORAM_ENABLE_VENDOR_TIMING 0
#endif

#if GIGADORAM_ENABLE_VENDOR_TIMING
#define DORAM_TIMING_ENABLED 1
#else
#define DORAM_TIMING_ENABLED 0
#endif

#if DORAM_TIMING_ENABLED
#define DORAM_TIMING_ADD(var, delta) \
    do {                              \
        (var) += (delta);             \
    } while (0)
#define DORAM_TIMING_SET(var, value) \
    do {                             \
        (var) = (value);             \
    } while (0)
#else
#define DORAM_TIMING_ADD(var, delta) \
    do {                              \
        (void)sizeof(var);            \
        (void)sizeof(delta);          \
    } while (0)
#define DORAM_TIMING_SET(var, value) \
    do {                             \
        (void)sizeof(var);           \
        (void)sizeof(value);         \
    } while (0)
#endif

