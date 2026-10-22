# Public probe configuration example

`probe.template.json` contains an invented demonstration probe, not a production
probe design. Its sequences and thresholds are illustrative and are not validated
for experimental use.

The web workflow and web AI share `<configured data_dir>/probes/*.yaml`.
Each filename is the method name (for example, `MiP-seq.yaml`), and each file
contains that method's configuration directly, without an outer method-name key.
The existing `/workflow/builtin_probes` endpoint still returns a method-to-config
mapping. Files are reread on every lookup; bad YAML files are logged and skipped
individually. The packaged CLI templates are independent and unchanged.

For a one-time conversion of an existing JSON library:

```sh
python -m uprobe.http.utils.migrate_probes /path/to/probe.json /path/to/data/probes
```

The conversion preserves all fields and refuses to overwrite an existing
directory. The original JSON remains available for backup, but is not read at
runtime. Customize the YAML files in the new directory after migration.
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
