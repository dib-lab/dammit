# Advanced Configuration

Dammit's overall memory and CPU usage can be specified at the command line.
The [annotation pipelines](pipelines.md) section contains info on the
recommended minimum resources for each pipeline.


## **`dammit config`**

```
Usage: dammit config [OPTIONS] COMMAND [ARGS]...

  Show dammit configuration information.

Options:
  --help  Show this message and exit.

Commands:
  busco-groups      Lists the available BUSCO group databases.
  clean-temp        Clear out shared dammit temp files.
  show-default      Show the selected default configuration file.
  show-directories  List dammit directory locations.
```



## Tool-Specific Specification

Tool-specific parameters can be modified via a custom configuration file.


## Distributing dammit Jobs across a Cluster







