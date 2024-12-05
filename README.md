# README


## Illumina NGS Aligner Benchmarking

This repository contains scripts and .Rmd/quarto files for benchmarking simulated reads from the E. coli K12 MG1655 genome using NCBI's ART read simulator. The primary script is called `benchmark_aligners.sh` and usage statement is available `-?|--help`. The benchmarking script is used best when the source is read and understood. In many cases, default parameters are used, and an important parameter to consider is `-p|--parallel` for parallelization. 


If you want to look at GNU time results for `exit` code and `real` seconds elapsed, check out the `${aligner}_timed.tsv` file for details on the distribution of real seconds elapsed using default parameters. Note that by using default parameters, the aligners may perform differently in terms of elapsed time, but this is not the primary output KPI we're looking at with aligners anyways.


The primary KPI are statistics like alignment rates, mate orientation, and proper pairing, which will in the future be diagnosed with `samtools` and `Picard`. For the time being, bare with me as the project develops.



# What is the rest of this garbage?

This is a README template I use for most projects, and the rest of the project has not been completed so the format will be useful in the future, should the project continue.

![Project Logo](logo_url.png)

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Build Status](https://travis-ci.org/username/project-name.svg?branch=master)](https://travis-ci.org/username/project-name)
[![Coverage Status](https://coveralls.io/repos/github/username/project-name/badge.svg?branch=master)](https://coveralls.io/github/username/project-name?branch=master)

A brief description of what this project does and who it's for.

## Table of Contents

- [Features](#features)
- [Demo](#demo)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [API Reference](#api-reference)
- [Contributing](#contributing)
- [Tests](#tests)
- [FAQ](#faq)
- [License](#license)
- [Contact](#contact)
- [Acknowledgements](#acknowledgements)

## Features

- Feature 1
- Feature 2
- Feature 3

## Demo

![Demo](demo.gif)

## Installation and Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
crate_name = "0.1.0"
```

Then, add this to your crate root (`main.rs` or `lib.rs`):

```rust
use crate_name;

fn main() {
    crate_name::some_function();
}
```

For more detailed examples, please refer to the [documentation](https://docs.rs/crate_name).
```


## Configuration

Explain how to configure your project, if applicable.

```json
{
  "key": "value",
  "anotherKey": "anotherValue"
}
```

## API Reference

### `functionName(param1, param2)`

Description of the function.

- `param1` (Type): Description of param1
- `param2` (Type): Description of param2

Returns: Description of return value

## Contributing

Contributions are always welcome! Please read the [contribution guidelines](CONTRIBUTING.md) first.

## Tests

```bash
npm test
```

## FAQ

**Q: Frequently asked question?**

A: Answer to the question.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

Matt Ralston - [Website](https://matthewralston.github.io/) - professional.bio.coder@gmail.com

Project Link: [https://github.com/YourUsername/template_py](https://github.com/username/template_py)

## Acknowledgements

- [Library/Resource 1](https://example.com)
- [Library/Resource 2](https://example.com)
- [Library/Resource 3](https://example.com)



<!--
Thanks of course to my fans (and haters). Yeah i see you.... but i dont.
Thanks to my former mentors Andrew S, Charles T, Brian C, Mark R, Isaac N, Carlos R, and my newer bosses Punita J and Kyle L.
Thanks to the Pap lab for the first dataset that I continue to use.
Thank you to Ryan for the food and stuff. I actually made this project specifically so you and I could converse...
Thanks to Blahah for tolerating someone snooping and imitating his Ruby style.
Thanks to Erin for getting my feet wet in this new field. You are my mvp.
Thanks to Rachel for the good memories and friendship. And Sophie too. veggies n' R love.
Thanks to Yasmeen for the usual banter.
Thanks to A for the newer banter.
Thanks to Max, Robin, and Robert for the good memories in St. Louis. What's new?
Thanks to Fred for the good memories. Hope you're on soon.
Thanks to Nichole for the cookies and good memories. And your cute furballs too! Hope you're well
Thanks to S for the lessons, convos, and even embarassing moments. You're kind of awesome to me.
Thanks to a few friends I met in 2023 that reminded me I have a lot to learn about friendship, dating, and street smarts.
Thanks to them even more now that I got it xd up.

Thanks to the people of NCC for the Doordash money. It might not be much but I don't have it twisted (I do.)



Thanks to D from BCCS.
Thanks to C4H&H. I'm 'healthier' now, but I really still think I need more support than just BCCS. it's urgent.
Thanks to CT and family. Your love and support means the world to me.
Thanks to AS and family. Your support made a difference. Praying for better employment and opportunities.

And thanks to my family and friends.
Go Blue Hens
-->
