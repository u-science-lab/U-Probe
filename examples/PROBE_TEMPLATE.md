# Public probe configuration example

`probe.template.json` contains an invented demonstration probe, not a production
probe design. Its sequences and thresholds are illustrative and are not validated
for experimental use.

Copy it to the configured data directory as `probe.json` and customize it locally.
Keep the production file outside the repository when it contains unpublished designs.

- Barcode lengths are optional. Add `length` to a barcode or `default_length` to
  `barcodes` when a fixed length is required.
- Attributes without a filter condition are calculated without implicit filtering.
- GC conditions use fractions from 0 to 1; the frontend displays percentages.
- Filter conditions support one-sided bounds and strict comparisons. Complex
  expressions remain unchanged and can be edited in the advanced expression field.
- Sort fields refer to the original attribute names and accept arrays of fields.

This file demonstrates configuration syntax only. It contains no production
barcode sequences or private probe structures.
