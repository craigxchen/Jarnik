# Known four-point family: Walsh diagnostic

## Purpose

This note checks the affine-cube/Walsh framework against the first known
four-point endpoint example. The example is not a counterexample to Walsh
rigidity. Instead, it exhibits exactly the permitted structure: an almost
character-compatible two-dimensional cube with three heavy character blocks and
one low-weight defect layer.

## Circle and points

The known example has

```text
N = 567454025
  = 5^2 * 61 * 233 * 1597.
```

Choose Gaussian factors

```text
5    = (1+2i)(1-2i),
61   = (5+6i)(5-6i),
233  = (8+13i)(8-13i),
1597 = (21+34i)(21-34i).
```

One endpoint arc of angular width `N^(-1/4)` contains the four points

```text
(23200,5405),
(23189,5452),
(23176,5507),
(23171,5528).
```

Using exponent coordinates in the displayed prime order, and absorbing the
necessary Gaussian units, their exponent vectors are

```text
A = (1,0,0,0),
B = (0,1,0,1),
C = (2,0,1,1),
D = (0,1,1,0).
```

## Binary threshold code

Expand the exponent of `5^2` into two threshold layers. In the layer order

```text
(5,1), (5,2), (61,1), (233,1), (1597,1),
```

the four binary words are

```text
A = 10000,
B = 00101,
C = 11011,
D = 00110.
```

The second `5`-layer has pattern

```text
0010.
```

Remove that one layer. The remaining words are

```text
A = 1000,
B = 0101,
C = 1011,
D = 0110.
```

Their bitwise XOR is zero, and after translating by `A` they form the affine
plane

```text
{0, u, v, u+v} subset F_2^4.
```

Thus the example is exactly a two-dimensional Boolean affine cube after one
low-weight defect layer is removed.

## Character classes

On the four cube vertices, the remaining layers restrict to the three nonzero
Walsh characters:

- the first `5`-layer and the `61`-layer are opposite orientations of the same
  character;
- the `233`-layer is the second character;
- the `1597`-layer is the third character.

The corresponding conductor weights are

```text
log(5*61),
log 233,
log 1597,
```

while the discarded defect has weight

```text
log 5.
```

Relative to

```text
H = log N,
```

the three character blocks carry approximately

```text
28%, 27%, 37%
```

of the conductor height, and the defect carries approximately `8%`.

Therefore all three active character blocks are heavy. This is exactly why the
four-point family survives the Walsh argument: the arithmetic-rigidity theorem
forces active character blocks to be heavy, but three heavy blocks can still
distinguish four cube vertices.

## Lesson

The diagnostic supports the current contradiction program:

1. global endpoint near-equality polarizes prime layers;
2. high additive energy produces Walsh-type affine cubes;
3. Gaussian arithmetic forces every active Walsh block to be heavy;
4. dimension five is impossible because `32` vertices require at least five
   active characters, but four heavy blocks are the maximum;
5. the known four-point example lies below that threshold and has precisely the
   expected heavy-block structure.
