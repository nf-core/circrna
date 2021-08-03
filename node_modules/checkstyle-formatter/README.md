# checkstyle-formatter

Simple Checkstyle data formatter. Formats data to an XML string, leaving the
reporting details to the user.

[![npm Version][npm-badge]][npm]
[![Build Status][build-badge]][build-status]
[![Test Coverage][coverage-badge]][coverage-result]
[![Dependency Status][dep-badge]][dep-status]

## Installation

Install using npm:

    $ npm install checkstyle-formatter

## Usage

```js
var checkstyleFormatter = require('checkstyle-formatter');
var results = [
    {
        filename: 'foo.js',
        messages: [
            {
                line: 1,
                column: 2,
                severity: 'warning',
                message: 'the quick'
            },
            {
                line: 3,
                column: 4,
                severity: 'error',
                message: 'brown fox'
            }
        ]
    },
    {
        filename: 'bar.js',
        messages: [
            {
                line: 5,
                column: 6,
                severity: 'warning',
                message: 'jumped over'
            },
            {
                line: 7,
                column: 8,
                severity: 'error',
                message: 'the lazy dog'
            }
        ]
    }
];

console.log(checkstyleFormatter(results));
// <?xml version="1.0" encoding="utf-8"?>
// <checkstyle version="4.3">
// <file name="foo.js">
// <error line="1" column="2" severity="warning" message="the quick" />
// <error line="3" column="4" severity="error" message="brown fox" />
// </file>
// <file name="bar.js">
// <error line="5" column="6" severity="warning" message="jumped over" />
// <error line="7" column="8" severity="error" message="the lazy dog" />
// </file>
// </checkstyle>
```

## Changelog

#### 1.1.0
- Add support for optional `source` information

#### 1.0.0
- Initial release

## License

MIT

[build-badge]: https://img.shields.io/travis/jimf/checkstyle-formatter/master.svg
[build-status]: https://travis-ci.org/jimf/checkstyle-formatter
[npm-badge]: https://img.shields.io/npm/v/checkstyle-formatter.svg
[npm]: https://www.npmjs.org/package/checkstyle-formatter
[coverage-badge]: https://img.shields.io/coveralls/jimf/checkstyle-formatter.svg
[coverage-result]: https://coveralls.io/r/jimf/checkstyle-formatter
[dep-badge]: https://img.shields.io/david/jimf/checkstyle-formatter.svg
[dep-status]: https://david-dm.org/jimf/checkstyle-formatter
