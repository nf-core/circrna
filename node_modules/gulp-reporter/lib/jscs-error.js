'use strict';
const LintError = require('./lint-error');
/**
 * JSCS错误对象
 *
 * @class JSCSError
 * @extends {LintError}
 */
class JSCSError extends LintError {
	/**
	 * Creates an instance of JSCSError.
	 *
	 * @param {Object} error JSCS的原始error对象
	 *
	 * @memberOf JSCSError
	 */
	constructor (error) {
		super({

			// 文件名
			fileName: error.filename,

			// JSCS无警告，错误等级全部算错误
			severity: 'error',

			// 行号
			lineNumber: error.line,

			// 列号
			columnNumber: error.column,

			// 错误信息
			// message: error.message,

			// 错误ID
			// rule: error.rule,

			// 源代码上下文
			source: error.element._sourceCode,

			// 报错插件
			plugin: 'JSCS',

			// 文档
			doc: error.rule && `http://jscs.info/rule/${error.rule}`,
		}, error, {
			// 错误信息
			message: error.message.replace(new RegExp('^' + error.rule + '\\:\\s*'), ''),
		});
	}
}

module.exports = JSCSError;
