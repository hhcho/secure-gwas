#include "crypto.h"
#include "errors.h"

#include <assert.h>
#include <openssl/sha.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int setup_prf_key(AES_KEY* key, const byte* buf, uint32_t buflen) {
  if (buflen != PRF_KEY_BYTES) {
    return ERROR_PRF_KEYLEN_INVALID;
  }

  AES_128_Key_Expansion(buf, key);

  return ERROR_NONE;
}

int generate_prf_key(AES_KEY* key) {
  byte keybuf[PRF_KEY_BYTES];

  FILE* f = fopen("/dev/urandom", "r");
  if (f == NULL) {
    return ERROR_RANDOMNESS;
  }

  int bytes_read = fread(keybuf, 1, PRF_KEY_BYTES, f);
  if (bytes_read != PRF_KEY_BYTES) {
    return ERROR_RANDOMNESS;
  }

  fclose(f);

  AES_128_Key_Expansion(keybuf, key);
  memset(keybuf, 0, sizeof(keybuf));

  return ERROR_NONE;
}
