'use strict';
const LintError = require('./lint-error');

const JSHINT_SEVERITY_MAP = {
	E: 'error',
	W: 'warn',
	I: 'info',
};

/**
 * JSHint错误对象
 *
 * @class JSHintError
 * @extends {LintError}
 */
class JSHintError extends LintError {
	/**
	 * Creates an instance of JSHintError.
	 *
	 * @param {any} error JSHint的原始error对象
	 *
	 * @memberOf JSHintError
	 */
	constructor (error) {
		super({
			// 文件名
			fileName: error.file,

			// JSHint错误等级
			severity: error.code && /^[EWI]/.test(error.code) ? JSHINT_SEVERITY_MAP[error.code[0]] : JSHINT_SEVERITY_MAP.E,

			// 行号
			lineNumber: error.line,

			// 列号
			columnNumber: error.character,

			// 错误信息
			message: error.reason,

			// 错误ID
			rule: error.code,

			// 源代码上下文
			source: error.evidence,

			// 报错插件
			plugin: 'JSHint',

		}, error);
	}
}

module.exports = JSHintError;
