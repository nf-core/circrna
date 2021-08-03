'use strict';
const BufferStreams = require('bufferstreams');
const cssBrowserReporter = require('./css-browser-reporter');
const jsBrowserReporter = require('./js-browser-reporter');
/**
 * 在文件内容上追加错误汇报信息
 * @param  {Vinyl}  file   需要检查的文件
 * @param  {Array}  errors 错误信息
 * @param  {Buffer} [buf=buffile.contents]  文件内容buffer, 默认从file.contents获取
 * @return {Buffer}        新的文件内容
 */
function appendReporter (file, errors, buf) {
	buf = buf || file.contents;

	// 在buffer中的代码中注入报错语句
	let contentReporter;

	if (file.postcss || file.csslint) {
		contentReporter = cssBrowserReporter(errors);
	} else if (file.jshint || file.eslint || file.jscs) {
		contentReporter = `\n;\n(${jsBrowserReporter})(${JSON.stringify(errors, circular())});`;
	}

	if (contentReporter) {
		return Buffer.concat([buf, Buffer.from(contentReporter)]);
	} else {
		return buf;
	}
}

/**
 * 为Vinyl对象追加错误汇报内容
 * @param  {Vinyl}  file     要处理的文件
 * @param  {Array}  [errors] 错误信息
 * @return {undefined}
 */
function browserReporter (file, errors) {
	errors = errors || file.report.errors;
	if (file.isStream()) {
		file.contents = file.contents.pipe(new BufferStreams((err, buf, done) => {
			if (err) {
				done(err, buf);
			} else {
				done(null, appendReporter(file, errors, buf));
			}
		}));
	} else if (file.isBuffer()) {
		file.contents = appendReporter(file, errors);
	}
}

function circular () {
	const cache = [];
	return (key, value) => {
		if (typeof value === 'object' && value !== null) {
			if (cache.indexOf(value) !== -1) {
				// Circular reference found, discard key
				return;
			}
			// Store value in our collection
			cache.push(value);
		}
		return value;
	};
}

module.exports = browserReporter;
