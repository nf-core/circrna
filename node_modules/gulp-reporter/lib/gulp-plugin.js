'use strict';
const through = require('through2');
const PluginError = require('plugin-error');
const getErrors = require('./get-errors');
const reporter = require('./reporter');
const ciReporter = require('./ci-reporter');
const getOptions = require('./get-options');

module.exports = function (options) {
	let files = [];

	options = getOptions(options);

	function transform (file, encoding, done) {
		const report = file.report || (file.report = {});
		options(file).then(options => {
			report.options = options;
			return getErrors(file, options).then(errors => {
				report.errors = errors;
				report.fail = reporter(file, options);
				files.push(file);
				done(null, file);
			});
		}).catch(error => {
			done(new PluginError('gulp-reporter', error), file);
		});
	}

	function flush (done) {
		const fails = files.filter(file => (
			file.report.options.fail && file.report.fail
		));

		ciReporter(files).then(() => {
			if (fails.length) {
				throw new PluginError('gulp-reporter', {
					message: 'Lint failed for: ' + fails.map(file => (
						file.relative.replace(/\\/g, '/')
					)).join(', '),
					showStack: false,
					showProperties: false,
				});
			}
		}).then(done, done);
		files = [];
	}

	return through.obj(transform, flush);
};
