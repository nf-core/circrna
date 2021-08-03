'use strict';
var fs = require('fs');
var path = require('path');
var gulpIgnore = require('gulp-ignore');

var appendStars = function (dirname, str) {
  var prefix = '';
  if (dirname) {
    prefix = dirname + '/';
  }
  return [
    prefix + str + '**',
    prefix + str + '/**'
  ];
};

module.exports = function (gitignorePath) {
  gitignorePath = path.resolve(gitignorePath || '.gitignore');
  var dirname = path.dirname(path.relative(process.cwd(), gitignorePath));
  if (dirname === '.') {
    dirname = '';
  }

  var contents = fs.readFileSync(gitignorePath, 'utf8');
  var ignoredFiles = contents.split('\n')
    .map(function (str) {
      return str.trim();
    })
    .filter(Boolean) // ignore empty lines
    .map(appendStars.bind(null, dirname))
    .reduce(function (m, paths) {
      return m.concat(paths);
    }, []);

  return gulpIgnore.exclude(ignoredFiles, {dot: true});
};
