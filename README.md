# ohash
Optimal hash string matching algorithms improve hash string matching algorithms.
The idea of the new algorithm is to consider substrings of length $q$
 when the pattern has no two substrings of length $q$ hashed to same the value.

It consists of three algorithms:
- ohash1.c: perfect hashing for $q=1$, hashing with possible collisions for $2\le q \le 10$ and HASH8 when $q>10$; 
- ohash2.c: perfect hashing for $q=1$ and $q=2$, hashing with possible collisions for $3\le q \le 10$ and HASH8 when $q>10$;
- ohash3.c: perfect hashing for $q=1$ and $q=2$, HASH3 for $3\le q \le 10$ and HASH8 when $q>10$.

The algorithms have been implemented so that they can directly be plugged in the
String Matching Algorithm Research Tool.
