# to-time
[![Build Status](https://travis-ci.org/hafuta/to-time.svg?branch=master)](https://travis-ci.org/hafuta/to-time) [![npm version](https://badge.fury.io/js/to-time.svg)](https://badge.fury.io/js/to-time)

Utility for converting textual time periods to time units (milliseconds, seconds, minutes, hours, etc..).

######Due to the lack of precision in floating point numbers arithmetic and the need of keeping the results precise the utility is using the [bignumber.js](https://github.com/MikeMcl/bignumber.js/) library for arithmetic operations on numbers.

## Install

First install the package and save it to package.json using npm:
```sh
npm install --save to-time
```

To require in the browser:

Both files located in _lib_ directory already include the bignumber.js dependency and there is no need to include this module in the browser.
```html
<!-- access using window.toTime -->
<script src="node_modules/to-time/lib/to-time.min.js"></script>
```

To require when using NodeJS:
```node
const toTime = require('to-time');
```

## Usage
Converting from textual time period to time units

```javascript
toTime('1 hour').seconds(); //3600
//same as:
toTime('1h').seconds(); //3600

toTime('1 Year 365 Days 4 Hours').hours(); //17524
//same as:
toTime('1y 365d 4h').hours(); //17524
```

Useful for usage in methods such as setInterval and setTimeout which consume the second argument in milliseconds, it is much clearer for someone who will read the code.
```javascript
//Instead of using 43200000 milliseconds (equivalent to 12 hours) we can do the following
setInterval(() => {
  //Do something here
}, toTime('12h').ms());

//Instead of using 5400000 milliseconds (equivalent to 1.5 hour)
setTimeout(() => {
  //Do something here
}, toTime.fromHours(1.5).ms());
```

##### Allowed suffixes (all case-insensetive)

* Year, Years, Y
* Week, Weeks, W
* Day, Days, D
* Hour, Hours, H
* Minute, Minutes, M
* Second, Seconds, S
* Millisecond, Milliseconds, MS


#### Initializing using factory methods
It is also possible to create a a TimeFrame instance by calling the static factory methods.
The result will be a TimeFrame object similar to the one that is created by invoking the function with a textual time period.
```javascript
toTime.fromHours(4).addMinutes(30).hours(); //4.5
toTime.fromYears(4).addWeeks(4).days(); //1488
```
###### Available factory methods
* fromMilliseconds
* fromSeconds
* fromMinutes
* fromHours
* fromDays
* fromWeeks
* fromYears


#### Appenders methods
It is possible to add additional units to the TimeFrame object by invoking one of the appender methods on the returned instance:

```javascript
toTime('0.5 hour').addMinutes(30).seconds(); //3600
toTime.fromHours(2).addMinutes(30).minutes(); //150
```

###### Available appenders:
* addMilliseconds
* addSeconds
* addMinutes
* addHours
* addDays
* addWeeks
* addYears

#### Getters methods
The getter methods are used to get the value of the TimeFrame object in a specific time unit.
* milliseconds : Number
* ms (alias to milliseconds)
* minutes : Number
* hours : Number
* days : Number
* weeks : Number
* years : Number
* humanize : String

#### Converting the TimeFrame object to human readable format
It is also possible to use to-time in order to convert time units into human readable format.
Example:
```javascript
const frame = toTime.fromMilliseconds(500005050505005);
frame.humanize(); //15855 Years, 2 Weeks, 6 Days, 11 Hours, 48 Minutes, 25 Seconds, 5 Milliseconds
```

## Contributing

##### Running tests

* Make sure to write tests, run new & existing tests using:
  ```shell
  npm run test
  ```
  
* Check for source code & tests code styling by running eslint:
  ```shell
  npm run lint
  ```

* If tests are passing and eslint doesn't return any error -> Create pull request

## To Do
- [ ]  Add Karma for testing in browser environment (Currently testing server side using Mocha)


## License
MIT
