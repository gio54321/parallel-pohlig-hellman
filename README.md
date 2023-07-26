# Parallel Pohlig-Hellman

This is a parallel Pohlig-Hellman + baby-step giant-step implementation made for the final project of _Parallel and distributed systems, paradigms and models_ course at unipi.

It employs a few tricks to make computation fast and somewhat scalable, namely:

- all the standard implementation tricks of sequential bsgs,
- shared hash table exploiting atomic operations,
- implicit storage in the table of the resulting element `g^i`, avoiding allocating big nums and improving memory usage.
