'use strict';

const POSTCSS_SEVERITY_MAP = {
	'warning': 'warn',
};
const LintError = require('./lint-error');

/**
 * PostCSS错误对象
 *
 * @class PostCSSError
 * @extends {LintError}
 */
class PostCSSError extends LintError {
	/**
	 * Creates an instance of PostCSSError.
	 *
	 * @param {Object} error 原始的PostCSS错误对象
	 *
	 * @memberOf PostCSSError
	 */
	constructor (error) {
		super({

			// 文件名
			fileName: error.file || (error.input && error.input.file),

			// 错误等级默认error，后面会覆盖
			// severity: 'error',

			// 行号
			lineNumber: error.line,

			// 列号
			columnNumber: error.column,

			// 错误信息
			message: error.text.replace(new RegExp('\\s*\\(' + error.rule + '\\)$'), ''),

			// 错误ID
			// rule: error.rule,

			// 源代码上下文
			// source: (error.node && error.node.type && error.node.type !== 'root') ? String(error.node) : '',

			// 报错插件
			plugin: 'PostCSS',
			doc: error.rule && error.plugin === 'stylelint' && `https://stylelint.io/user-guide/rules/${error.rule}/`,
		}, error, {
			// 错误等级
			severity: (error.severity && (POSTCSS_SEVERITY_MAP[error.severity] || error.severity)) || 'error',
		});
	}
}

module.exports = PostCSSError;
