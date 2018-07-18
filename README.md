# XiP

Xoshiro in Python (XiP) is a library for fast generation of pseudorandom numbers in python.
It uses a C++ implementation of two xoshiro generators.
The library mirrors the interface of python's random module.
This means that projects using python stock randomness can easily be ported to this library.

XiP offers two xoshiro style generators: xoshiro256** and xoshiro256+.
By default it uses xoshiro256**, which is the best choice for most applications.
xoshiro256+ should only be used to generate floating point numbers.
See [here](http://xoshiro.di.unimi.it/) for additional information.

## Advantages

xoshiro generators offer faster generation of pseudorandom numbers than Mersenne Twister (python's default PRNG) while passing all BigCrush tests.

You can achieve significant speedup by bundling your pseudorandom number generation into a single call to XiP and then use these numbers one after another.

## Usage

Use the provided Makefile to compile XiP's C++ lib into a shared object.
The simply `import xip` in python to use the XiP class.

Creating more than one XiP object allows you to access multiple independent streams of pseudorandomness simultaneously.

**Example**
~~~
from xip import XiP

gen = XiP() # If no seed is provided, then XiP initializes using urandom
X = gen.random(100)

... # use X to feed your application with randomness
~~~

## License

See the LICENSE file for license rights and limitations (MIT).
