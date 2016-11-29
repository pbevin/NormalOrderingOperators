# NormalOrderingOperators

https://github.com/pbevin/NormalOrderingOperators

Pete Bevin <pete@petebevin.com>

## What is this?

It converts a product of quantum fields into *normal order*, meaning that
all the creation operators end up to the left of all the annihilation
operators.

Obviosuly, you probably aren't going to be interested in this unless
you do quantum field theory for fun :)

You are free to copy, modify, and distribute NormalOrderingOperators
with attribution under the terms of the MIT license. See the LICENSE.md
file for details.


## How to use it

Before installing NormalOrderingOperators, you need:

  * Git
  * Python 3

```
$ git clone https://github.com/pbevin/NormalOrderingOperators.git
$ cd NormalOrderingOperators
$ python3 NOO.py 'ad1s[k1].ad1s[k2].ah1s[p1mu].cd1s[p1m].ad1s[p2m].cd1s[q1].cd1s[q2].ch1s[q3]'
```

## Contributing

Contributions are welcome as bug reports, feature requests, or pull requests.

### Bugs and feature requests

Please use the [issue tracker](https://github.com/pbevin/NormalOrderingOperators/issues)
to report bugs or request features.

### Pull Requests

Follow the normal Github process for pull requests, and please include tests with
new features.

### Testing

Most functions are unit tested using `doctest`. You can run the suite with the `--selftest` flag:

```
$ python3 NOO.py --selftest
OK
```

You can also run `test.sh`, which also runs a regression test.
