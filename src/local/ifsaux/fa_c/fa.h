#pragma once

#include "lfi_type.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void faitou64_ (integer64 * KREP, integer64 * KNUMER, logical * LDNOMM, character * CDNOMF, character * CDSTTU,   
                       logical * LDERFA, logical * LDIMST, integer64 * KNIMES, integer64 * KNBARP, integer64 * KNBARI,      
                       character * CDNOMC, character_len CDNOMF_len, character_len CDSTTU_len, character_len CDNOMC_len);

extern void fairme64_ (integer64 * KREP, integer64 * KNUMER, character * CDSTTU, character_len CDSTTU_len);

extern void facies64_ (character * CDNOMC, integer64 * KTYPTR, real64 * PSLAPO, real64 * PCLOPO, real64 * PSLOPO,  
                       real64 * PCODIL, integer64 * KTRONC, integer64 * KNLATI, integer64 * KNXLON, integer64 * KNLOPA,       
                       integer64 * KNOZPA, real64 * PSINLA, integer64 * KNIVER, real64 * PREFER, real64 * PAHYBR,       
                       real64 * PBHYBR, logical *LDGARD, character_len CDNOMC_len);

#ifdef __cplusplus
}
#endif


