A library and example executables for generating developable surfaces attached
to a provided C²-smooth R³-space curve, where a parameterization allows
steering the surface normal across the length of the curve. No arc-length
parameterization is required for the space curve.

The resultant surface can be exported as a 3D mesh in `.obj` and as its
flattened 2D pattern in `.svg` format. Compiles to a standalone plain
`wasm32-wasip1` target for ~200kB of binary size, requiring nothing but
stdin/stdout communication.

```bash
cargo build --target wasm32-wasip1  -p dc-plot --release
```
