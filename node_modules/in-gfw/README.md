in-gfw
=====
[![NPM version](https://img.shields.io/npm/v/in-gfw.svg?style=flat-square)](https://www.npmjs.com/package/in-gfw)
[![Travis](https://img.shields.io/travis/gucong3000/in-gfw.svg?&label=Linux)](https://travis-ci.org/gucong3000/in-gfw)
[![AppVeyor](https://img.shields.io/appveyor/ci/gucong3000/in-gfw.svg?&label=Windows)](https://ci.appveyor.com/project/gucong3000/in-gfw)
[![Codecov](https://img.shields.io/codecov/c/github/gucong3000/in-gfw.svg)](https://codecov.io/gh/gucong3000/in-gfw)
[![David](https://img.shields.io/david/gucong3000/in-gfw.svg)](https://david-dm.org/gucong3000/in-gfw)

Identify if current location is located in mainland China.

## Install

```bash
npm install in-gfw
```

## Usage

```js
const inGFW = require("in-gfw");
inGFW().then(console.log);	// `true` for located in mainland China
inGFW.os().then(console.log);	// `true` for system located in mainland China
inGFW.net().then(console.log);	// `true` for network located in mainland China
```

## API

```js
inGFW(blockedHost, cnHost);
inGFW.sync(blockedHost, cnHost);	// Synchronous version of `inGFW()`
```
Get result by `inGFW.os()` and fallback to `inGFW.net()`

```js
inGFW.os();
inGFW.osSync();	// Synchronous version of `inGFW.os()`
```

- Windows: Check if current location settings is `PRC`.
  > Control Panel: Regional and language -> Location
- POSIX systems: Check if [timezone](https://en.wikibooks.org/wiki/Puredyne/Date_and_Timezone) is set to `Beijing`, `Chongqing`, `Shanghai`, `Urumqi` or `PRC`.

```js
inGFW.net(blockedHost, cnHost);
inGFW.netSync(blockedHost, cnHost);	// Synchronous version of `inGFW.net()`
```
Based on the speed of network access to identify if current location is located in mainland China.

- `blockedHost`

  Type: `string|URL`

  Default: `"www.npmjs.com"`

  host for speed test that blocked by GFW

- `cnHost`

  Type: `string|URL`

  Default: `"npm.taobao.org"`

  host for speed test that mirrored in mainland China.
