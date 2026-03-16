pub mod obj;
pub mod svg;

#[must_use = "This must be written to a file to do anything."]
pub struct StrFileData {
    pub contents: String,
}

/// This is for use in `std::io::Write as _`.
impl core::ops::Deref for StrFileData {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.contents.as_bytes()
    }
}

/// This is for use in `std::fs::write`.
impl AsRef<[u8]> for StrFileData {
    fn as_ref(&self) -> &[u8] {
        self.contents.as_bytes()
    }
}
