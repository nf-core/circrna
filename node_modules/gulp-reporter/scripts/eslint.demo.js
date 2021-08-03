'use strict';
const vfs = require('vinyl-fs');
const eslint = require('gulp-eslint');
const reporter = require('../');
const path = require('path');

vfs.src('test/fixtures/eslint/invalid.js', {
	cwd: path.join(__dirname, '..'),
})
	.pipe(eslint())
	.pipe(reporter({
		// author: '1',
		// blame: false,
	})).on('error', error => {
		process.exitCode = 1;
		if (error) {
			console.error(String(error));
		}
	});
