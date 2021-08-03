'use strict';

const LintError = require('./lint-error');
const CSSLINT_SEVERITY_MAP = {
	'warning': 'warn',
};

/**
 * CSSLint错误对象
 *
 * @class CSSLintError
 * @extends {LintError}
 */
class CSSLintError extends LintError {
	/**
	 * Creates an instance of CSSLintError.
	 *
	 * @param {Object} error 原始的CSSLint错误对象
	 *
	 * @memberOf CSSLintError
	 */
	constructor (error) {
		const rule = error.rule || {};
		super({
			// 文件名
			// fileName: error.fileName,

			// 错误等级
			severity: CSSLINT_SEVERITY_MAP[error.type] || error.type,

			// 行号
			lineNumber: error.line,

			// 列号
			columnNumber: error.col,

			// 源代码上下文
			source: error.evidence,

			// 报错插件
			plugin: 'CSSLint',

			// 文档
			doc: rule.url,
		}, error, {
			// 错误ID
			rule: rule.id,
		});
	}
}

module.exports = CSSLintError;
