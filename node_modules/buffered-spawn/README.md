# buffered-spawn

[![NPM version][npm-image]][npm-url] [![Downloads][downloads-image]][npm-url] [![Build Status][travis-image]][travis-url] [![Dependency status][david-dm-image]][david-dm-url] [![Dev Dependency status][david-dm-dev-image]][david-dm-dev-url]

[npm-url]:https://npmjs.org/package/buffered-spawn
[downloads-image]:http://img.shields.io/npm/dm/buffered-spawn.svg
[npm-image]:http://img.shields.io/npm/v/buffered-spawn.svg
[travis-url]:https://travis-ci.org/IndigoUnited/node-buffered-spawn
[travis-image]:http://img.shields.io/travis/IndigoUnited/node-buffered-spawn/master.svg
[david-dm-url]:https://david-dm.org/IndigoUnited/node-buffered-spawn
[david-dm-image]:https://img.shields.io/david/IndigoUnited/node-buffered-spawn.svg
[david-dm-dev-url]:https://david-dm.org/IndigoUnited/node-buffered-spawn#info=devDependencies
[david-dm-dev-image]:https://img.shields.io/david/dev/IndigoUnited/node-buffered-spawn.svg

Buffered child_process#spawn.


## Installation

`$ npm install buffered-spawn`


## Why

- Easy to use
- Uses [cross-spawn](http://github.com/IndigoUnited/node-cross-spawn) by default, which fixes windows [issues](https://github.com/joyent/node/issues/2318)
- Supports callback & promise style


## Usage

In terms of arguments, they are equal to node's [spawn](http://nodejs.org/api/child_process.html#child_process_child_process_spawn_command_args_options).

```js
const bufferedSpawn = require('buffered-spawn');

// Callback style
bufferedSpawn('git', ['clone', 'git@github.com/bower/bower'], { cwd: '.' }, (err, stdout, stderr) => {
    if (err) {
        // Both stdout and stderr are also set on the error object
        return console.error(`Command failed with error code of #${err.status}`);
    }

    console.log(stdout);
    console.log(stderr);
});

// ...or Promise style
bufferedSpawn('git', ['clone', 'git@github.com/bower/bower'], { cwd: '.' })
.then((output) => {
    console.log(output.stdout);
    console.log(output.stderr);
}, (err) => {
    // Both stdout and stderr are also set on the error object
    console.error(`Command failed with error code of #${err.status}`);
});
```

The actual child process is available if necessary:

```js
const buffspawn('buffered-spawn');

// Callback style
const cp = buffspawn('git', ['clone', 'git@github.com/bower/bower'], () => {}};

// ...or Promise style
const promise = buffspawn('git', ['clone', 'git@github.com/bower/bower']);
const cp = promise.cp;
```

As said before, `buffered-spawn` uses `cross-spawn` to actually spawn the process. If you are having trouble running Windows such as [wmic](https://msdn.microsoft.com/en-us/library/bb742610.aspx) which has its own arguments parser, you may disable like so:

```js
const cp = buffspawn('wmic', [
    'logicaldisk', 'where', 'DeviceID="Z:"',
    'get' 'freeSpace,size'
], { crossSpawn: false }, () => {}};
```

## Tests

`$ npm test`


## License

Released under the [MIT License](http://www.opensource.org/licenses/mit-license.php).
