use core::intrinsics;

#[cfg(feature = "mem")]
use mem::{memcpy, memmove, memset};

// NOTE This function and the ones below are implemented using assembly because they using a custom
// calling convention which can't be implemented using a normal Rust function
#[naked]
#[cfg_attr(not(test), no_mangle)]
pub unsafe fn __aeabi_uidivmod() {
    asm!("push {lr}
          sub sp, sp, #4
          mov r2, sp
          bl __udivmodsi4
          ldr r1, [sp]
          add sp, sp, #4
          pop {pc}");
    intrinsics::unreachable();
}

#[naked]
#[cfg_attr(not(test), no_mangle)]
pub unsafe fn __aeabi_uldivmod() {
    asm!("push {r4, lr}
          sub sp, sp, #16
          add r4, sp, #8
          str r4, [sp]
          bl __udivmoddi4
          ldr r2, [sp, #8]
          ldr r3, [sp, #12]
          add sp, sp, #16
          pop {r4, pc}");
    intrinsics::unreachable();
}

#[naked]
#[cfg_attr(not(test), no_mangle)]
pub unsafe fn __aeabi_idivmod() {
    asm!("push {r0, r1, r4, lr}
          bl __divsi3
          pop {r1, r2}
          muls r2, r2, r0
          subs r1, r1, r2
          pop {r4, pc}");
    intrinsics::unreachable();
}

#[naked]
#[cfg_attr(not(test), no_mangle)]
pub unsafe fn __aeabi_ldivmod() {
    asm!("push {r4, lr}
          sub sp, sp, #16
          add r4, sp, #8
          str r4, [sp]
          bl __divmoddi4
          ldr r2, [sp, #8]
          ldr r3, [sp, #12]
          add sp, sp, #16
          pop {r4, pc}");
    intrinsics::unreachable();
}

#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_dadd(a: f64, b: f64) -> f64 {
    ::float::add::__adddf3(a, b)
}

#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_fadd(a: f32, b: f32) -> f32 {
    ::float::add::__addsf3(a, b)
}

#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_dsub(a: f64, b: f64) -> f64 {
    ::float::sub::__subdf3(a, b)
}

#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_fsub(a: f32, b: f32) -> f32 {
    ::float::sub::__subsf3(a, b)
}

#[cfg(not(all(feature = "c", target_arch = "arm", not(target_os = "ios"), not(thumbv6m))))]
#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_idiv(a: i32, b: i32) -> i32 {
    ::int::sdiv::__divsi3(a, b)
}

#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_lasr(a: i64, b: u32) -> i64 {
    ::int::shift::__ashrdi3(a, b)
}

#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_llsl(a: u64, b: u32) -> u64 {
    ::int::shift::__ashldi3(a, b)
}

#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_llsr(a: u64, b: u32) -> u64 {
    ::int::shift::__lshrdi3(a, b)
}

#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_lmul(a: u64, b: u64) -> u64 {
    ::int::mul::__muldi3(a, b)
}

#[cfg(not(all(feature = "c", target_arch = "arm", not(target_os = "ios"), not(thumbv6m))))]
#[cfg_attr(not(test), no_mangle)]
pub extern "aapcs" fn __aeabi_uidiv(a: u32, b: u32) -> u32 {
    ::int::udiv::__udivsi3(a, b)
}

#[cfg(not(feature = "c"))]
#[cfg_attr(not(test), no_mangle)]
pub extern "C" fn __aeabi_ui2d(a: u32) -> f64 {
    ::float::conv::__floatunsidf(a)
}

// TODO: These aeabi_* functions should be defined as aliases
#[cfg(not(feature = "mem"))]
extern "C" {
    fn memcpy(dest: *mut u8, src: *const u8, n: usize) -> *mut u8;
    fn memmove(dest: *mut u8, src: *const u8, n: usize) -> *mut u8;
    fn memset(dest: *mut u8, c: i32, n: usize) -> *mut u8;
}

// FIXME: The `*4` and `*8` variants should be defined as aliases.

#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memcpy(dest: *mut u8, src: *const u8, n: usize) {
    memcpy(dest, src, n);
}
#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memcpy4(dest: *mut u8, src: *const u8, n: usize) {
    memcpy(dest, src, n);
}
#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memcpy8(dest: *mut u8, src: *const u8, n: usize) {
    memcpy(dest, src, n);
}

#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memmove(dest: *mut u8, src: *const u8, n: usize) {
    memmove(dest, src, n);
}
#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memmove4(dest: *mut u8, src: *const u8, n: usize) {
    memmove(dest, src, n);
}
#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memmove8(dest: *mut u8, src: *const u8, n: usize) {
    memmove(dest, src, n);
}

// Note the different argument order
#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memset(dest: *mut u8, n: usize, c: i32) {
    memset(dest, c, n);
}
#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memset4(dest: *mut u8, n: usize, c: i32) {
    memset(dest, c, n);
}
#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memset8(dest: *mut u8, n: usize, c: i32) {
    memset(dest, c, n);
}

#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memclr(dest: *mut u8, n: usize) {
    memset(dest, 0, n);
}
#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memclr4(dest: *mut u8, n: usize) {
    memset(dest, 0, n);
}
#[cfg_attr(not(test), no_mangle)]
pub unsafe extern "aapcs" fn __aeabi_memclr8(dest: *mut u8, n: usize) {
    memset(dest, 0, n);
}
