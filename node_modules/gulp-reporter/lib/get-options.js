'use strict';
const hasCheckstyle = require('./has-checkstyle');
const gitAuthor = require('./git-author');
const termSize = require('term-size');
const log = require('fancy-log');
const toTime = require('to-time');
const ci = require('ci-info');
function unixTimestamp (now) {
	return Math.floor((now || Date.now()) / 1000);
}

module.exports = function (options) {
	const authorCache = {};
	const termColumns = ci.isCI
		? 160
		: Math.max(termSize().columns, 80);

	function getAuthor (cwd) {
		return authorCache[cwd] || (authorCache[cwd] = gitAuthor(cwd));
	}

	function getOptions (file) {
		if (typeof options === 'function') {
			options = options(file);
		}
		options = Object.assign({
			maxLineLength: 512,
			browser: false,
			output: !(ci.isCI && (ci.APPVEYOR || ci.CIRCLE || (ci.JENKINS && hasCheckstyle))),
			blame: true,
			fail: true,
			sort: true,
		}, options);
		return Promise.resolve(options.blame && getAuthor(file.cwd)).then(author => {
			options = Object.assign({
				author,
			}, options);

			if (typeof options.author === 'string') {
				if (/@/.test(options.author)) {
					options.author = {
						email: options.author,
					};
				} else {
					options.author = {
						name: options.author,
					};
				}
			}

			let expires = options.expires;

			if (expires) {
				if (typeof expires === 'string') {
					try {
						expires = toTime(expires).seconds();
					} catch (ex) {
						expires = new Date(expires);
					}
				}
				if (typeof expires === 'number') {
					if (expires <= 0 || isNaN(expires)) {
						throw new TypeError('`options.expires` must be greater than 0.');
					}
					options._expiresTime = ((author && author.time) || unixTimestamp()) - expires;
				} else if (expires.getTime) {
					expires = expires.getTime();
					if (isNaN(expires)) {
						throw new TypeError('`options.expires` must be valid `Date`.');
					}
					options._expiresTime = unixTimestamp(expires);
				} else {
					throw new TypeError('`options.expires` must be `Number`, `Date` or `string`.');
				}
			}

			options._termColumns = termColumns;

			let writable = options.output;
			if (writable) {
				if (typeof writable.write === 'function') {
					writable = writable.write.bind(writable);
				} else if (typeof writable !== 'function') {
					writable = log.warn;
				}
				options.output = writable;
			}
			return options;
		});
	}
	return getOptions;
};
