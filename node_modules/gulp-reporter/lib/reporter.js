'use strict';
const browserReporter = require('./browser-reporter');
const formatter = require('./formatter');

function reporter (file, options) {
	function isFail () {
		if (typeof options.fail === 'function') {
			return error => options.fail(error, file);
		} else {
			return error => !error.demote && (!error.severity || error.severity === 'error');
		}
	}
	const errors = file.report.errors;
	if (!errors.length) {
		return;
	}
	const writable = options.output;
	if (writable) {
		writable(formatter(file, options));
	}
	if (options.browser) {
		browserReporter(file);
	}
	return errors.some(isFail());
}
module.exports = reporter;
