'use strict';
const demoteErrors = require('./demote-errors');
const sortError = require('./sort-errors');
const gitBlame = require('./git-blame');
const shortDocUrl = require('./short-doc-url');
const addPostcssSource = require('./add-postcss-source');
const get = require('lodash.get');

const linterMap = {
	standard: 'eslint',
};

/**
 * 提取`file`中的postcss、eslint、tslint、jscs、jshint等错误
 *
 * @param {Vinyl}  file    需要检查的文件
 * @param {Object} options 参数，这里要使用options.sort和options.filter
 * @returns {Promise} <Array>
 */
function getErrors (file, options) {
	// 开始收集错误
	let errors = [];
	function addErrors (path, hook) {
		let errs = get(file, path);
		if (!errs || ((hook && !(errs = hook(errs))) || !errs.length)) {
			return;
		}
		path = path.replace(/^(\w+).*$/, (s, linter) => {
			linter = linterMap[linter] || linter;
			return './' + linter + '-error';
		});
		const LintError = require(path);
		errs.forEach(error => errors.push(new LintError(error)));
	}

	addErrors('csslint.report.messages');
	addErrors('editorconfig.errors');
	addErrors('eslint.messages');
	addErrors('htmlhint.messages', messages => messages.map(message => message.error || message));
	addErrors('jscs.errors', errors => errors.getErrorList && errors.getErrorList());
	addErrors('jshint.results', results => results.map(result => result.error || result));
	addErrors('jsonlint.message', message => [message]);
	addErrors('postcss.messages', addPostcssSource.bind(this, file));
	addErrors('standard.results[0].messages');
	addErrors('stylint.msg', msg => msg.replace(/(\n*Stylint:\s+\d+\s+\w+\.)+$/i, '').split(/\n\n/g));
	addErrors('tslint.failures');

	// 所务收集完成
	if (!errors.length) {
		return Promise.resolve(errors);
	}

	// 获取git blame信息

	return (options.blame ? gitBlame(file).catch(() => {

	}) : Promise.resolve()).then(blames => {
		// blame信息写入error对象
		errors.forEach((error) => {
			if (!error.fileName) {
				error.fileName = file.path;
			}
			const blame = blames && error.lineNumber && blames[error.lineNumber];
			if (blame) {
				error.blame = blame.rev;
				error.source = blame.content;
			}
		});

		// 错误所在行的长度超出阈值时，删除错误
		const maxLineLength = options.maxLineLength;
		if (maxLineLength) {
			errors = errors.filter(error => {
				return !(error.columnNumber > maxLineLength || (error.source && /^.*$/m.exec(error.source)[0].length > maxLineLength));
			});
		}

		// 错误降级
		if (blames) {
			demoteErrors(errors, options);
		}

		// 错误排序
		if (options.sort) {
			errors = sortError(errors, typeof options.sort === 'function' ? options.sort : null);
		}

		errors = shortDocUrl(errors);

		if (options.mapper) {
			errors = errors.then(errors => (
				Promise.resolve(options.mapper(file)).then(mapper => (
					errors.map(mapper)
				)).then(errors => (
					Promise.all(errors)
				)).then(errors => (
					errors.filter(Boolean)
				))
			));
		}

		return errors;
	});
}

module.exports = getErrors;
