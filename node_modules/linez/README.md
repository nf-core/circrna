# linez

> Parses lines from text, preserving line numbers, offsets and line endings.

[![Build Status](https://secure.travis-ci.org/jedmao/linez.svg?branch=master)](http://travis-ci.org/jedmao/linez)
[![Dependency Status](https://gemnasium.com/jedmao/linez.svg)](https://gemnasium.com/jedmao/linez)
[![npm version](https://badge.fury.io/js/linez.svg)](http://badge.fury.io/js/linez)
[![codecov](https://codecov.io/gh/jedmao/linez/branch/master/graph/badge.svg)](https://codecov.io/gh/jedmao/linez)

[![npm](https://nodei.co/npm/linez.png?downloads=true)](https://nodei.co/npm/linez/)


## Getting Started

### Installation

```bash
$ npm install linez
```

### Importing

#### TypeScript

```ts
import * as linez from 'linez';
```

#### Babel/ES6+

```js
import linez from 'linez';
```

#### JavaScript

```js
var linez = require('linez');
```

## Introduction

By default, linez uses `/\r?\n/g` as the regular expression to detect newline character sequences and split lines. This regular expression is tuned for performance and only covers the most common newline types (i.e., `\n` and `\r\n`). If you have need for more newline character sequences, you can configure linez with the `configure` method.

```js
linez.configure({
  newlines: ['\n', '\r\n', '\r', '\u000B']
});
```

Setting this property will automatically create a piped regular expression for you and use it in any future `linez()` calls. You can make up your own newlines if you want. Linez doesn't care one way or the other.

```js
linez.configure({
  newlines: ['foo', 'bar']
});
```

This would be converted into `/(foo|bar)/g`. Newlines are just strings. They can be anything. There are, however, some known newline character sequences. Should you need them, refer to the following table:

| String   | Unicode        | Name                        |
| -------- |:-------------- |:--------------------------- |
| `\n`     | U+000A         | Line feed                   |
| `\r\n`   | U+000D, U+000A | Carriage Return + Line Feed |
| `\r`     | U+000D         | Carriage Return             |
| `\u000B` | U+000B         | Vertical Tab                |
| `\u000C` | U+000C         | Form Feed                   |
| `\u0085` | U+0085         | Next Line                   |
| `\u2028` | U+2028         | Line Separator              |
| `\u2029` | U+2029         | Paragraph Separator         |


### Byte Order Marks

Also referred to as BOM signatures, these are the bytes at the beginning of a file that indicating the encoding in which the file is written. Currently, linez only reads BOMs to detect the encoding and does not take into account the contents of the file.

#### Supported BOMs

- utf-8-bom
- utf-16le
- utf-16be

#### Unsupported BOMs

- utf-32le
- utf-32be

If linez detects an unsupported BOM, an error will be thrown, indicating that decoding the detected charset is not supported.

#### Default decoding

By default, the document will attempt to be decoded as utf8. This is the default behavior of [the Node API's conversion from buffers into strings](https://nodejs.org/api/buffer.html#buffer_buf_tostring_encoding_start_end).


# API


### configure(options: IOptions)

Configures linez to use the supplied options. Currently, only the newlines property is available, where you can specify any number of newline character sequences.

```js
linez.configure({
  newlines: ['\n', '\r\n', '\r', '\u000B']
});
```

### resetConfiguration()

Resets the configuration to the default settings, using `/\r?\n/g` as the newlines regular expression.


### Document

```ts
constructor(public lines: Line[]);
```

Calling the `toString()` method converts the document's lines into a string, discarding information about line numbers and offsets.


### Line

```ts
interface Line {
  offset: number;
  number: number;
  text: string;
  ending: string;
}
```

### Options

```ts
interface Options {
  newlines?: string[];
}
```


### linez(file: string|Buffer): Document

Parses text into a `Document`.

[The specs](https://github.com/jedmao/linez/blob/master/lib/linez.spec.ts) show some great usage examples.

```ts
var lines = linez('foo\nbar\nbaz').lines;
lines[1].offset; // 4
lines[1].number; // 2
lines[1].text; // bar
lines[1].ending; // \n
```

Note: You can also pass-in a Buffer.


## License

Released under the MIT license.
