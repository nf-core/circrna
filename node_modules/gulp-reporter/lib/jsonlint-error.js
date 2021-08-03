'use strict';
const LintError = require('./lint-error');
/**
 * JSON Lint错误对象
 *
 * @class JsonLintError
 * @extends {LintError}
 */
class JsonLintError extends LintError {
	/**
	 * Creates an instance of JsonLintError.
	 *
	 * @param {string} message JSON Lint的原始error消息
	 *
	 * @memberOf JsonLintError
	 */
	constructor (message) {
		message = message.match(/^Parse error on line (\d+):\n(.*)\n(-*\^)\n(.+)$/i);
		super({

			// JSON Lint无警告，错误等级全部算错误
			severity: 'error',

			// 行号
			lineNumber: +message[1],

			// 列号
			columnNumber: message[3].length,

			// 源代码上下文
			source: message[2],

			// 错误信息
			message: message[4],

			// 报错插件
			plugin: 'JSONLint',
		});
	}
}

module.exports = JsonLintError;
