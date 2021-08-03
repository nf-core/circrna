'use strict';
const HTMLHINT_SEVERITY_MAP = {
	'warning': 'warn',
};

const LintError = require('./lint-error');
/**
 * HTMLHint错误对象
 *
 * @class HTMLHintError
 * @extends {LintError}
 */
class HTMLHintError extends LintError {
	/**
	 * Creates an instance of HTMLHintError.
	 *
	 * @param {Object} error eslint原生错误对象
	 *
	 * @memberOf HTMLHintError
	 */
	constructor (error) {
		super({
			// 错误等级
			severity: HTMLHINT_SEVERITY_MAP[error.type] || error.type,

			// 行号
			lineNumber: error.line,

			// 列号
			columnNumber: error.col,

			// 错误信息
			// message: error.message,

			// 触发错误的规则
			// rule: error.rule.id,

			// 源代码上下文
			source: error.raw,

			// 报错插件
			plugin: 'HTMLHint',

			// 文档
			doc: error.rule.link,
		}, error, {
			// 触发错误的规则
			rule: error.rule.id,
		});
	}
}

module.exports = HTMLHintError;
